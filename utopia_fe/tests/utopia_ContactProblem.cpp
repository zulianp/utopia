#include "utopia_ContactProblem.hpp"
#include "utopia_assemble_contact.hpp"
#include "utopia_fe_core.hpp"
#include "utopia_Socket.hpp"

using namespace std;
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
	void assemble_elasticity(
		U &u, 
		LibMeshFEContext<SystemType> &context,
		DSMatrixd &stiffness_matrix,
		DVectord &force
		)
	{
		const int dim = u.size();

		double mu = 10., lambda = 10.;
		// auto e  = transpose(grad(u)) + grad(u); //0.5 moved below -> (2 * 0.5 * 0.5 = 0.5)
		// auto b_form = integral((mu * 0.5) * dot(e, e) + lambda * dot(div(u), div(u)));
		auto b_form = integral(dot(grad(u), grad(u)));

		DenseVector<Real> vec(dim);
		vec.zero();

		auto f 	    = vec_coeff(vec);
		auto l_form = integral(dot(f, u));

		auto ass = make_assembly([&]() -> void {
			double t = MPI_Wtime();

			assemble(u, u, b_form, l_form, *context.system.matrix, *context.system.rhs);

			t = MPI_Wtime() - t;

			printf("--------------------------------\n");
			printf("Assembly: %g seconds\n", t);
			printf("--------------------------------\n");
		});

		context.system.attach_assemble_object(ass);
		context.equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 1;
		context.equation_systems.solve();

		convert( *context.system.matrix, stiffness_matrix);
		convert( *context.system.rhs, force);

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
	}

	void ContactProblem::init_material()
	{
		if(mesh->mesh_dimension() == 2) {
			init_material_2d();
		} else {
			init_material_3d();
		}
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

		DSMatrixd D_inv = diag(d_inv);
		transfer_operator = D_inv * coupling;
		transfer_operator += local_identity(local_size(d).get(0), local_size(d).get(0));
		gap = D_inv * weighted_gap;

		if(comm.is_alone()) plot_scaled_normal_field(*mesh, normals, gap);
	}

	void ContactProblem::step()
	{
		disp("n_dofs:");
		disp(size(force));
		compute_contact_conditions();

		DVectord sol_c = local_zeros(local_size(force));
		DVectord rhs_c = transpose(orthogonal_trafo) * transpose(transfer_operator) * force;
		DSMatrixd K_c  = transpose(orthogonal_trafo) *
							DSMatrixd(
								transpose(transfer_operator) * 
								stiffness_matrix * 
								transfer_operator) * 
							orthogonal_trafo;


		std::shared_ptr< LinearSolver<DSMatrixd, DVectord> > linear_solver;

		// auto s_sol = size(sol_c);
		// if(s_sol.get(0) > 2.5e5) {
		// 	auto cg = std::make_shared<ConjugateGradient<DSMatrixd, DVectord> >();
		// 	cg->atol(1e-12);
		// 	cg->rtol(1e-12);
		// 	cg->stol(1e-12);
		// 	linear_solver = cg;
		// } else {
			linear_solver = std::make_shared<Factorization<DSMatrixd, DVectord> >();
		// }

		SemismoothNewton<DSMatrixd, DVectord> newton(linear_solver);
		newton.verbose(true);
		newton.set_active_set_tol(1e-8);
		newton.max_it(40);

		newton.set_box_constraints(make_upper_bound_constraints(make_ref(gap)));
		newton.solve(K_c, rhs_c, sol_c);

		displacement = transfer_operator * (orthogonal_trafo * sol_c);
		apply_displacement(displacement);
	}

	void ContactProblem::apply_displacement(const DVectord &displacement)
	{
		//FIXME
		// int sys_num = 0;

		// displaced_mesh = std::shared_ptr<MeshBase>(mesh->clone().release());
		// Read<DVectord> r_d(displacement);

		// for(auto n_it = displaced_mesh->local_nodes_begin(); n_it != displaced_mesh->local_nodes_end(); ++n_it) {
		// 	for(unsigned int c = 0; c < displaced_mesh->mesh_dimension(); ++c) {
		// 		const int dof_id = (*n_it)->dof_number(sys_num, c, 0);
		// 		(**n_it)(c) += displacement.get(dof_id);
		// 	}
		// }
	}

	void ContactProblem::save(const std::string &output_dir)
	{
		convert(displacement, *context_ptr->system.solution);
		ExodusII_IO(*context_ptr->mesh).write_equation_systems (output_dir + "/sol.e", context_ptr->equation_systems);

		convert(is_contact_node, *context_ptr->system.solution);
		ExodusII_IO(*context_ptr->mesh).write_equation_systems (output_dir + "/is_c_n.e", context_ptr->equation_systems);
	}
}
