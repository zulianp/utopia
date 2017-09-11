#include "utopia_ContactProblem.hpp"
#include "utopia_assemble_contact.hpp"
#include "utopia_fe_core.hpp"
#include "utopia_Socket.hpp"
#include "libmesh/parameter_vector.h"
#include "libmesh/parameter_accessor.h"

#include <memory>

//using namespace std;
using std::make_shared;
using std::shared_ptr;

using namespace libMesh;

namespace utopia {

	void ContactProblem::init(
		const libMesh::LibMeshInit &init, 
		const std::shared_ptr<Mesh> &mesh,
		const std::shared_ptr< ElasticityBoundaryConditions > &bc_ptr,
		std::vector< std::pair<int, int> > contact_pair_tags,
		double search_radius
		)
	{
		this->bc_ptr = bc_ptr;
		this->contact_pair_tags = contact_pair_tags;
		this->search_radius = search_radius;
		this->mesh = mesh;
		comm.set_mpi_comm(init.comm().get());
		init_discretization();
		init_material();
		iteration = 0;
		linear_solver = std::make_shared<Factorization<DSMatrixd, DVectord> >();
	}

	void ContactProblem::init_discretization()
	{
		context_ptr = make_shared<FEContextT>(mesh);

		spaces.clear();
		for(int i = 0; i < mesh->mesh_dimension(); ++i) {
			spaces.push_back( make_shared<FESpaceT>( fe_space(LAGRANGE, Order(1), *context_ptr) ) );
		}
	}

	template<class U, class SystemType>
	void assemble_mass_matrix(
		U &u, 
		LibMeshFEContext<SystemType> &context,
		DSMatrixd &mass_matrix)
	{
		const int dim = u.size();

		double mu = 10., lambda = 10.;
		// auto e  = transpose(grad(u)) + grad(u); //0.5 moved below -> (2 * 0.5 * 0.5 = 0.5)
		// auto b_form = integral((mu * 0.5) * dot(e, e) + lambda * dot(div(u), div(u)));
		auto b_form = integral(dot(u, u));
		
		long n_local_dofs = u.get(0).dof_map().n_local_dofs();
		mass_matrix = local_sparse(n_local_dofs, n_local_dofs, 20);
		assemble(u, u, b_form, mass_matrix, false);
	}

	template<class U, class SystemType>
	void assemble_elasticity(
		U &u, 
		LibMeshFEContext<SystemType> &context,
		DSMatrixd &stiffness_matrix,
		DVectord &force
		)
	{
		const int dim = u.size();

		double mu = 0.2, lambda = 0.2;
		auto e  = transpose(grad(u)) + grad(u); //0.5 moved below -> (2 * 0.5 * 0.5 = 0.5)
		auto b_form = integral((mu * 0.5) * dot(e, e) + lambda * dot(div(u), div(u)));
		// auto b_form = integral(dot(grad(u), grad(u)));

		DenseVector<Real> vec(dim);
		vec.zero();

		auto f 	    = vec_coeff(vec);
		auto l_form = integral(dot(f, u));

		auto ass = make_assembly([&]() -> void {
			double t = MPI_Wtime();

			assemble(u, u, b_form, l_form, *context.system.matrix, *context.system.rhs, false);

			t = MPI_Wtime() - t;

			printf("--------------------------------\n");
			printf("Assembly: %g seconds\n", t);
			printf("--------------------------------\n");
		});



		context.system.attach_assemble_object(ass);
		context.equation_systems.parameters.template set<unsigned int>("linear solver maximum iterations") = 1;
		context.equation_systems.solve();

		convert( *context.system.matrix, stiffness_matrix);
		convert( *context.system.rhs, force);
		apply_boundary_conditions(u.get(0), stiffness_matrix, force);
	}

	void ContactProblem::init_material_2d()
	{
		const int dim = 2;
		std::cout << "assmebling 2d problem" << std::endl;

		auto ux = fe_function(*spaces[0]);
		auto uy = fe_function(*spaces[1]);
		auto u = prod(ux, uy);

		bc_ptr->apply(ux, uy);

		context_ptr->equation_systems.init();
		ux.set_quad_rule(make_shared<libMesh::QGauss>(dim, SECOND));
		uy.set_quad_rule(make_shared<libMesh::QGauss>(dim, SECOND));
		
		assemble_elasticity(u, *context_ptr, stiffness_matrix, force);
		assemble_mass_matrix(u, *context_ptr, mass_matrix);
	}

	void ContactProblem::init_material_3d()
	{
		const int dim = 3;
		std::cout << "assmebling 3d problem" << std::endl;

		auto ux = fe_function(*spaces[0]);
		auto uy = fe_function(*spaces[1]);
		auto uz = fe_function(*spaces[2]);
		auto u = prod(ux, uy, uz);

		bc_ptr->apply(ux, uy, uz);

		context_ptr->equation_systems.init();
		ux.set_quad_rule(make_shared<libMesh::QGauss>(dim, SECOND));
		uy.set_quad_rule(make_shared<libMesh::QGauss>(dim, SECOND));
		uz.set_quad_rule(make_shared<libMesh::QGauss>(dim, SECOND));

		assemble_elasticity(u, *context_ptr, stiffness_matrix, force);
		assemble_mass_matrix(u, *context_ptr, mass_matrix);
	}

	void ContactProblem::init_material()
	{
		if(mesh->mesh_dimension() == 2) {
			init_material_2d();
		} else {
			init_material_3d();
		}

		old_displacement = local_zeros(local_size(force));

		velocity = local_zeros(local_size(force));
		old_velocity = local_zeros(local_size(force));

		acceleration = local_zeros(local_size(force));
		old_acceleration = local_zeros(local_size(force));

		internal_force = local_zeros(local_size(force));
		total_displacement =  local_zeros(local_size(force));
		
	}

	void ContactProblem::compute_contact_conditions()
	{

		unsigned int variable_number = 0;

		assemble_contact(
			comm,
			mesh, 
			utopia::make_ref(context_ptr->system.get_dof_map()), 
			variable_number, 
			coupling, 
			orthogonal_trafo, 
			weighted_gap, 
			normals,
			is_contact_node, 
			search_radius,
			contact_pair_tags,
			true);

		DVectord d = sum(coupling, 1);
		DVectord d_inv = local_zeros(local_size(d));

		{
			Write<DVectord> w_(d_inv);

			each_read(d, [&d_inv](const SizeType i, const double value) {
				if(value < -1e-8) {
					std::cerr << "negative el for " << i << std::endl;
				}

				if(std::abs(value) > 1e-15) {
					d_inv.set(i, 1./value);
				} else {
					d_inv.set(i, 1.);
				}
			});
		}

		boundary_mass_inv = diag(d_inv);
		transfer_operator = boundary_mass_inv * coupling;
		transfer_operator += local_identity(local_size(d).get(0), local_size(d).get(0));
		gap = boundary_mass_inv * weighted_gap;

		if(comm.is_alone()) plot_scaled_normal_field(*mesh, normals, gap, "time_series_r/r" + std::to_string(iteration));
	}

	void ContactProblem::step(const double dt)
	{
		disp("n_dofs:");
		disp(size(force));
		const int dim = mesh->mesh_dimension();

		// DVectord inertia = mass_matrix * acceleration;
		// apply_zero_boundary_conditions(spaces[0]->dof_map(), inertia); 

		DVectord rhs = (force - dt * internal_force);// - (dt*dt) * inertia;

		if(iteration > 0)
		{
			search_radius = 1e-8;
			auto r = range(old_displacement);
			for(auto i = r.begin(); i < r.end(); i += dim) {
				double lenght = 0.0;

				for(int j = 0; j < dim; ++j) {
					auto v =  old_displacement.get(i + j);
					lenght += v*v;
				}

				search_radius = std::max(search_radius, std::sqrt(lenght));
			}
		}

		search_radius += 1e-8;
		comm.all_reduce(&search_radius, 1, moonolith::MPIMax());
		disp("predicted search radius");
		disp(search_radius);

		compute_contact_conditions();

		DVectord sol_c = local_zeros(local_size(force));
		DVectord rhs_c = transpose(orthogonal_trafo) * transpose(transfer_operator) * rhs;
		DSMatrixd K_c  = transpose(orthogonal_trafo) *
							DSMatrixd(
								transpose(transfer_operator) * 
								stiffness_matrix * 
								transfer_operator) * 
							orthogonal_trafo;

		//predictor step
		//DVectord displacement_pred = displacement + dt * velocity;
		//DVectord displacement_pred_c = transpose(orthogonal_trafo) * displacement_pred
		//DVectord addmissible_pred = min(displacement_pred_c, gap)
		//DVectord pred = addmissible_pred - displacement_pred_c;

		SemismoothNewton<DSMatrixd, DVectord> newton(linear_solver);
		newton.set_active_set_tol(1e-8);
		newton.max_it(40);

		newton.set_box_constraints(make_upper_bound_constraints(make_ref(gap)));
		newton.solve(K_c, rhs_c, sol_c);

		displacement = transfer_operator * (orthogonal_trafo * sol_c);
		
		total_displacement += displacement;
		internal_force = stiffness_matrix * total_displacement;

		apply_zero_boundary_conditions(spaces[0]->dof_map(), internal_force);
		
		normal_stress = local_zeros(local_size(displacement));

		DVectord unscaled_stress = boundary_mass_inv * (force - dt * internal_force);
		{
			Write<DVectord> w_ns(normal_stress);
			Read<DVectord> r_n(normals);
			Read<DVectord> r_d(unscaled_stress);


			auto r = range(unscaled_stress);

			
			for(auto i = r.begin(); i < r.end(); i += dim) {
				double ns = 0.;
				for(int j = 0; j < dim; ++j) {
					ns += normals.get(i + j) * unscaled_stress.get(i + j);
					// ns += unscaled_stress.get(i + j) * unscaled_stress.get(i + j);
				}

				normal_stress.set(i, ns);
			}
		}


		apply_displacement(displacement);

		old_acceleration = acceleration;
		old_velocity = velocity;

		velocity = displacement - old_displacement;
		velocity *= 1./dt;

		acceleration = velocity - old_velocity;
		acceleration *= 1./dt;

		old_velocity = velocity;
		old_displacement = displacement;

		++iteration;
	}

	void ContactProblem::apply_displacement(const DVectord &displacement)
	{
		//FIXME
		int sys_num = 0;
		Read<DVectord> r_d(displacement);

		auto m_it  = mesh->local_nodes_begin();
		auto m_end = mesh->local_nodes_end();


		for(; m_it != m_end; ++m_it) { //, ++d_it)
			for(unsigned int c = 0; c < mesh->mesh_dimension(); ++c) {
				const int dof_id = (*m_it)->dof_number(sys_num, c, 0);
				(**m_it)(c) += displacement.get(dof_id);
			}
		}

		const int dim = mesh->mesh_dimension();
		std::vector<double> normal_stress_x(local_size(normal_stress).get(0)/dim);
		
		{
			auto r = range(normal_stress);
			Read<DVectord> r_n(normal_stress);
			
			for(auto i = r.begin(); i < r.end(); i += dim) {
				normal_stress_x[i/dim] = normal_stress.get(i);
			}
		}

		if(comm.is_alone()) plot_mesh_f(*mesh, &normal_stress_x[0], "time_series_m/m" + std::to_string(iteration));
	}

	void ContactProblem::save(const std::string &output_dir)
	{
		convert(total_displacement, *context_ptr->system.solution);
		ExodusII_IO(*mesh).write_equation_systems (output_dir + "/sol_" + std::to_string(iteration) + ".e", context_ptr->equation_systems);

		convert(is_contact_node, *context_ptr->system.solution);
		ExodusII_IO(*mesh).write_equation_systems (output_dir + "/is_c_n_" + std::to_string(iteration) + ".e", context_ptr->equation_systems);

		convert(normal_stress, *context_ptr->system.solution);
		ExodusII_IO(*mesh).write_equation_systems (output_dir + "/ns_" + std::to_string(iteration) + ".e", context_ptr->equation_systems);

		convert(internal_force, *context_ptr->system.solution);
		ExodusII_IO(*mesh).write_equation_systems (output_dir + "/if_" + std::to_string(iteration) + ".e", context_ptr->equation_systems);
	
		convert(acceleration, *context_ptr->system.solution);
		ExodusII_IO(*mesh).write_equation_systems (output_dir + "/a_" + std::to_string(iteration) + ".e", context_ptr->equation_systems);

	}
}
