#include "utopia_ContactProblem.hpp"
#include "utopia_assemble_contact.hpp"
#include "utopia_fe_core.hpp"
#include "utopia_Socket.hpp"
#include "libmesh/parameter_vector.h"
#include "libmesh/parameter_accessor.h"
#include "libmesh/nemesis_io.h"
#include "moonolith_profiler.hpp"
#include "utopia_NormalTangentialCoordinateSystem.hpp"

#include <memory>

using std::make_shared;
using std::shared_ptr;

using namespace libMesh;

namespace utopia {

	ContactProblem::ContactProblem()
	{
		verbose = true;
		dynamic_contact = false;
		is_inpulse_ = false;
	}

	void ContactProblem::init(
		const libMesh::LibMeshInit &init, 
		const std::shared_ptr<MeshBase> &mesh,
		const std::shared_ptr< ElasticityBoundaryConditions > &bc_ptr,
		const std::shared_ptr< ElasticityForcingFunction > &ff_ptr,
		std::vector< std::pair<int, int> > contact_pair_tags,
		double search_radius
		)
	{
		if(verbose) std::cout << "Contact problem: initializing" << std::endl;

		MOONOLITH_EVENT_BEGIN("init");

		this->bc_ptr = bc_ptr;
		this->ff_ptr = ff_ptr;
		this->contact_pair_tags = contact_pair_tags;
		this->search_radius = search_radius;
		this->mesh = mesh;
		comm.set_mpi_comm(init.comm().get());
		init_discretization();
		init_material();
		iteration = 0;
		linear_solver = std::make_shared<Factorization<DSMatrixd, DVectord> >();
		// linear_solver = std::make_shared<BiCGStab<DSMatrixd, DVectord> >();
		// linear_solver = std::make_shared<ConjugateGradient<DSMatrixd, DVectord, HOMEMADE> >();
		output = std::make_shared<Nemesis_IO>(*mesh);
		// output = std::make_shared<ExodusII_IO>(*mesh);
		MOONOLITH_EVENT_END("init");
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
		DSMatrixd &mass_matrix,
		DVectord &inverse_mass_vector)
	{
		std::cout << "Contact problem: assembling mass matrix..." << std::flush;
		
		const int dim = u.size();
		auto b_form = integral(dot(u, u));
		
		long n_local_dofs = u.get(0).dof_map().n_local_dofs();

		double t = 0.0;
		auto ass = make_assembly([&]() -> void {
			t = MPI_Wtime();
			assemble_no_constraints(u, u, b_form, *context.system.matrix);
			t = MPI_Wtime() - t;
		});

		context.system.attach_assemble_object(ass);
		context.equation_systems.parameters.template set<unsigned int>("linear solver maximum iterations") = 1;
		context.equation_systems.solve();

		convert(*context.system.matrix, mass_matrix);
		DVectord lumped = sum(mass_matrix, 1);
		mass_matrix = diag(lumped);

		inverse_mass_vector = 1./lumped;

		std::cout << "done: (" << t << " seconds)" << std::endl; 
	}

	template<class U, class SystemType>
	void assemble_elasticity(
		U &u, 
		LibMeshFEContext<SystemType> &context,
		ContactProblem::ElasticityForcingFunction &ff,
		DSMatrixd &stiffness_matrix,
		DVectord &external_force
		)
	{
		std::cout << "Contact problem: assembling stiffness matrix..." << std::flush;
		
		const int dim = u.size();

		// auto mu     = block_var(6., {{1, 10.}, {4, 5.}});
		// auto lambda = block_var(6., {{1, 10.}, {4, 5.}});

		double mu = 730;
		double lambda = 376;
		// double mu = 10;
		// double lambda = 10;

		auto e  = transpose(grad(u)) + grad(u); //0.5 moved below -> (2 * 0.5 * 0.5 = 0.5)
		auto b_form = integral((0.5 * mu) * dot(e, e) + lambda * dot(div(u), div(u)));

		DenseVector<Real> vec_0(dim);
		vec_0.zero();

		DenseVector<Real> vec_1(dim);
		ff.fill(vec_1);

		auto f_0    = vec_coeff(vec_0);
		auto f_1    = vec_coeff(vec_1);
		auto l_form = integral(dot(f_0, u)) + integral(dot(f_1, u), ff.block_id());

		double t = 0.;
		auto ass = make_assembly([&]() -> void {
			t = MPI_Wtime();
			assemble(u, u, b_form, l_form, *context.system.matrix, *context.system.rhs, false);
			t = MPI_Wtime() - t;
		});

		context.system.attach_assemble_object(ass);
		context.equation_systems.parameters.template set<unsigned int>("linear solver maximum iterations") = 1;
		context.equation_systems.solve();

		convert(*context.system.matrix, stiffness_matrix);
		convert(*context.system.rhs, external_force);

		apply_boundary_conditions(u.get(0), stiffness_matrix, external_force);

		std::cout << "done: (" << t << " seconds)" << std::endl; 
		// write("s.m", stiffness_matrix);
	}

	void ContactProblem::init_material_2d()
	{
		const int dim = 2;
		if(verbose) std::cout << "assmebling 2d problem" << std::endl;

		auto ux = fe_function(*spaces[0]);
		auto uy = fe_function(*spaces[1]);
		auto u = prod(ux, uy);

		bc_ptr->apply(ux, uy);

		context_ptr->equation_systems.init();
		ux.set_quad_rule(make_shared<libMesh::QGauss>(dim, SECOND));
		uy.set_quad_rule(make_shared<libMesh::QGauss>(dim, SECOND));
		
		assemble_elasticity(u,  *context_ptr, *ff_ptr, stiffness_matrix, external_force);
		assemble_mass_matrix(u, *context_ptr, mass_matrix, inverse_mass_vector);
	}

	void ContactProblem::init_material_3d()
	{
		const int dim = 3;
		if(verbose) std::cout << "assmebling 3d problem" << std::endl;

		auto ux = fe_function(*spaces[0]);
		auto uy = fe_function(*spaces[1]);
		auto uz = fe_function(*spaces[2]);
		auto u = prod(ux, uy, uz);

		bc_ptr->apply(ux, uy, uz);

		context_ptr->equation_systems.init();
		ux.set_quad_rule(make_shared<libMesh::QGauss>(dim, SECOND));
		uy.set_quad_rule(make_shared<libMesh::QGauss>(dim, SECOND));
		uz.set_quad_rule(make_shared<libMesh::QGauss>(dim, SECOND));

		assemble_elasticity(u,  *context_ptr, *ff_ptr, stiffness_matrix, external_force);
		assemble_mass_matrix(u, *context_ptr, mass_matrix, inverse_mass_vector);
	}

	void ContactProblem::init_aux_system()
	{
		//init aux system for plotting
		auto &aux = context_ptr->equation_systems.add_system<libMesh::LinearImplicitSystem>("aux");
		
		var_num_aux.push_back( aux.add_variable("disp_x", Order(1), LAGRANGE) );
		var_num_aux.push_back( aux.add_variable("disp_y", Order(1), LAGRANGE) );

		const int dim = mesh->mesh_dimension();

		if(dim > 2) 
			var_num_aux.push_back( aux.add_variable("disp_z", Order(1), LAGRANGE) );

		var_num_aux.push_back( aux.add_variable("vel_x", Order(1), LAGRANGE) );
		var_num_aux.push_back( aux.add_variable("vel_y", Order(1), LAGRANGE) );

		if(dim > 2) 
			var_num_aux.push_back( aux.add_variable("vel_z", Order(1), LAGRANGE) );

		var_num_aux.push_back( aux.add_variable("f_x", Order(1), LAGRANGE) );
		var_num_aux.push_back( aux.add_variable("f_y", Order(1), LAGRANGE) );

		if(dim > 2) 
			var_num_aux.push_back( aux.add_variable("f_z", Order(1), LAGRANGE) );


		var_num_aux.push_back( aux.add_variable("c_x", Order(1), LAGRANGE) );
		var_num_aux.push_back( aux.add_variable("c_y", Order(1), LAGRANGE) );

		if(dim > 2) 
			var_num_aux.push_back( aux.add_variable("c_z", Order(1), LAGRANGE) );

		var_num_aux.push_back( aux.add_variable("is_contact_boundary", Order(1), LAGRANGE) );

		var_num_aux.push_back( aux.add_variable("rays_x", Order(1), LAGRANGE) );
		var_num_aux.push_back( aux.add_variable("rays_y", Order(1), LAGRANGE) );

		if(dim > 2) 
			var_num_aux.push_back( aux.add_variable("rays_z", Order(1), LAGRANGE) );

		aux.init();	
		aux.update();

	}

	void ContactProblem::assemble_velocities()
	{
		Write<DVectord> w(velocity);

		libMesh::DenseVector<libMesh::Real> v(mesh->mesh_dimension());

		std::vector<libMesh::dof_id_type> indices;
		for(auto e_it = mesh->active_local_elements_begin(); 
			e_it != mesh->active_local_elements_end(); 
			++e_it) {

			spaces[0]->dof_map().dof_indices(*e_it, indices);
			v.resize(indices.size());
			v.zero();

			iv_ptr->fill_velocity((*e_it)->subdomain_id(), v);
			
			for(uint i = 0; i < indices.size(); ++i) {
				velocity.set(indices[i], v(i));
			}
		}
	}

	void ContactProblem::init_material()
	{
		const int dim = mesh->mesh_dimension();
		if(dim == 2) {
			init_material_2d();
		} else {
			init_material_3d();
		}

		old_displacement_increment = local_zeros(local_size(external_force));
		displacement_increment =  local_zeros(local_size(external_force));

		velocity = local_zeros(local_size(external_force));

		if(iv_ptr) {
			assemble_velocities();
		}

		internal_force = local_zeros(local_size(external_force));
		total_displacement =  local_zeros(local_size(external_force));	

		DVectord selector = local_values(local_size(external_force).get(0), 1.);
		apply_zero_boundary_conditions(spaces[0]->dof_map(), selector);
		internal_mass_matrix = diag(selector) * mass_matrix;
		constrained_mass_matrix = mass_matrix;
		set_identity_at_constraint_rows(spaces[0]->dof_map(), constrained_mass_matrix);

		is_contact_node =  local_zeros(local_size(external_force));
		rays = local_zeros(local_size(external_force));
		normals = local_zeros(local_size(external_force));
		gap = local_zeros(local_size(external_force));

		init_aux_system();
	}

	void ContactProblem::compute_contact_conditions()
	{
		if(verbose) {
			comm.barrier();
			std::cout << "Contact problem: compute_contact_conditions" << std::endl;
		}

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

		if(comm.is_alone() && utopia::Utopia::Instance().get("plot") == "true") {
			plot_scaled_normal_field(*mesh, normals, gap, "time_series_r/r" + std::to_string(iteration));
		}
	}


	void ContactProblem::implicity_euler(const double dt)
	{
		DSMatrixd &K = stiffness_matrix;
		DVectord rhs = dt * (external_force - internal_force);

		DVectord sol_c = local_zeros(local_size(external_force));
		DVectord rhs_c = transpose(orthogonal_trafo) * transpose(transfer_operator) * rhs;
		DSMatrixd K_c  = transpose(orthogonal_trafo) *
		DSMatrixd(
			transpose(transfer_operator) * 
			K * 
			transfer_operator) * 
		orthogonal_trafo;


		// SemismoothNewton<DSMatrixd, DVectord, PETSC_EXPERIMENTAL> newton(linear_solver);
		SemismoothNewton<DSMatrixd, DVectord> newton(linear_solver);
		newton.verbose(true);
		newton.max_it(40);

		newton.set_box_constraints(make_upper_bound_constraints(make_ref(gap)));
		newton.solve(K_c, rhs_c, sol_c);

		displacement_increment = transfer_operator * (orthogonal_trafo * sol_c);		
		total_displacement += displacement_increment;

		new_internal_force = stiffness_matrix * total_displacement;
		apply_zero_boundary_conditions(spaces[0]->dof_map(), new_internal_force);
	}

	void ContactProblem::classic_newmark(const double dt)
	{
		DVectord &u_old = total_displacement;
		DVectord u_older = u_old - old_displacement_increment;
		DVectord rhs = 4. * ( 1./(dt*dt) * ( internal_mass_matrix * (u_old - u_older) ) + external_force - stiffness_matrix * (3./4. * u_old + 1./4. * u_older) );
		DSMatrixd K  = 4./(dt*dt) * internal_mass_matrix + stiffness_matrix;

		SemismoothNewton<DSMatrixd, DVectord> newton(linear_solver);
		newton.verbose(true);
		newton.max_it(40);

		DVectord dummy = local_values(local_size(rhs).get(0), 1000000);
		newton.set_box_constraints(make_upper_bound_constraints(make_ref(dummy)));

		displacement_increment = local_zeros(local_size(rhs));
		newton.solve(K, rhs, displacement_increment);

		total_displacement += displacement_increment;
		new_internal_force = stiffness_matrix * total_displacement;
		apply_zero_boundary_conditions(spaces[0]->dof_map(), new_internal_force);

		velocity = (1/dt) * displacement_increment;
	}

	void ContactProblem::classic_newmark_with_contact(const double dt)
	{
		DVectord &u_old = total_displacement;
		DVectord u_older = u_old - old_displacement_increment;
		DVectord rhs = 4. * ( 1./(dt) * ( internal_mass_matrix * velocity ) + external_force - stiffness_matrix * (3./4. * u_old + 1./4. * u_older) );
		DSMatrixd K  = 4./(dt*dt) * internal_mass_matrix + stiffness_matrix;

		DVectord sol_c = local_zeros(local_size(external_force));
		DVectord rhs_c = transpose(orthogonal_trafo) * transpose(transfer_operator) * rhs;
		DSMatrixd K_c  = transpose(orthogonal_trafo) *
		DSMatrixd(
			transpose(transfer_operator) * 
			K * 
			transfer_operator) * 
		orthogonal_trafo;

		SemismoothNewton<DSMatrixd, DVectord> newton(linear_solver);
		newton.verbose(true);
		newton.max_it(40);

		newton.set_box_constraints(make_upper_bound_constraints(make_ref(gap)));
		newton.solve(K_c, rhs_c, sol_c);

		displacement_increment = transfer_operator * (orthogonal_trafo * sol_c);	

		total_displacement += displacement_increment;
		new_internal_force = stiffness_matrix * total_displacement;
		apply_zero_boundary_conditions(spaces[0]->dof_map(), new_internal_force);

		velocity = (1./dt) * displacement_increment;

		//compute acceleration
		// DVectord dt_x_M_x_acc = (dt/2.) * (2. * external_force - internal_force - new_internal_force);
		// DVectord vel_inc = local_zeros(local_size(dt_x_M_x_acc));
		// apply_zero_boundary_conditions(spaces[0]->dof_map(), dt_x_M_x_acc);
		// // solve(mass_matrix, dt_x_M_x_acc, vel_inc);
		// vel_inc = e_mul(inverse_mass_vector, dt_x_M_x_acc);

		// velocity += vel_inc;
	}


	void ContactProblem::classic_newmark_with_contact_2(const double dt)
	{
		DVectord &u_old = total_displacement;
		DVectord u_older = u_old - old_displacement_increment;
		DVectord pred = orthogonal_trafo * utopia::min(orthogonal_trafo * (dt * velocity), gap);
		
		DVectord rhs = 4. * ( 1./(dt*dt) * ( internal_mass_matrix * pred ) + external_force - stiffness_matrix * (3./4. * u_old + 1./4. * u_older) );
		DSMatrixd K  = 4./(dt*dt) * internal_mass_matrix + stiffness_matrix;

		DVectord sol_c = local_zeros(local_size(external_force));
		DVectord rhs_c = transpose(orthogonal_trafo) * transpose(transfer_operator) * rhs;
		DSMatrixd K_c  = transpose(orthogonal_trafo) *
		DSMatrixd(
			transpose(transfer_operator) * 
			K * 
			transfer_operator) * 
		orthogonal_trafo;

		SemismoothNewton<DSMatrixd, DVectord> newton(linear_solver);
		newton.verbose(true);
		newton.max_it(40);

		newton.set_box_constraints(make_upper_bound_constraints(make_ref(gap)));
		newton.solve(K_c, rhs_c, sol_c);

		displacement_increment = transfer_operator * (orthogonal_trafo * sol_c);	

		total_displacement += displacement_increment;
		new_internal_force = stiffness_matrix * total_displacement;
		apply_zero_boundary_conditions(spaces[0]->dof_map(), new_internal_force);

		velocity = (1./dt) * displacement_increment;


		// 2.* ((_disp[_qp] - _disp_old[_qp]) / (_dt * _dt) - _vel_old[_qp] / _dt -
		//                         accel_old * (0.5 - _beta));


	}

	void ContactProblem::classic_newmark_beta(const double dt)
	{

	}

	static void make_block_system(DSMatrixd &K_11, DSMatrixd &K_21, DSMatrixd &K_22, DSMatrixd &K,
								  DVectord &rhs_1, DVectord &rhs_2, DVectord &rhs)
	{
		Size s = local_size(K_11);
	 	K = local_sparse(s.get(0) * 2, s.get(1) * 2, 60);
		Range rr = row_range(K_11);

		rhs = local_zeros(s.get(0) * 2);

		{
			Read<DSMatrixd> r_11(K_11), r_21(K_21), r_22(K_22);
			Write<DSMatrixd> w_K(K);

			for(auto i = rr.begin(); i < rr.end(); ++i) {
				RowView<DSMatrixd> rw_11(K_11, i);
				for(auto j = 0; j < rw_11.n_values(); ++j) {
					K.set(i, rw_11.get_col_at(j),  rw_11.get_value_at(j));
				}

				RowView<DSMatrixd> rw_21(K_21, i);
				for(auto j = 0; j < rw_21.n_values(); ++j) {
					K.set(i, rr.extent() + rw_21.get_col_at(j),  rw_21.get_value_at(j));
				}

				RowView<DSMatrixd> rw_22(K_22, i);
				for(auto j = 0; j < rw_22.n_values(); ++j) {
					K.set(rr.extent() + i, rr.extent() + rw_22.get_col_at(j), rw_22.get_value_at(j));
				}
			}
		}

		{
			Read<DVectord> r_1(rhs_1), r_2(rhs_2);
			Write<DVectord> w_K(rhs);

			for(auto i = rr.begin(); i < rr.end(); ++i) {
				rhs.set(i, rhs_1.get(i));
				rhs.set(rr.extent() + i, rhs_2.get(i));
			}
		}
	}

	static void unmake_block_system(DVectord &sol, DVectord &sol_1, DVectord &sol_2)
	{
		auto rr = range(sol_1);

		Read<DVectord> r_K(sol);
		Write<DVectord> w_1(sol_1), w_2(sol_2);
		
		for(auto i = rr.begin(); i < rr.end(); ++i) {
			sol_1.set(i, sol.get(i));
			sol_2.set(i, sol.get(rr.extent() + i));
		}
	}

	void ContactProblem::contact_stabilized_newmark_monolithic(const double dt)
	{
		DVectord displacement_increment_pred = dt * velocity;
		DVectord displacement_increment_pred_nt = utopia::min(orthogonal_trafo * displacement_increment_pred, gap);
		displacement_increment_pred = orthogonal_trafo * displacement_increment_pred_nt;
		

		DVectord rhs_1 = mass_matrix * displacement_increment_pred + (dt*dt/2.) * (external_force - internal_force);
		DSMatrixd K_11 = internal_mass_matrix + (dt*dt/4.) * stiffness_matrix;
		DSMatrixd K_21 = (dt*dt/2.) * stiffness_matrix;
		DSMatrixd &K_22 = constrained_mass_matrix;

		DVectord rhs_2 = internal_mass_matrix * velocity + dt/2. * (external_force - internal_force);
		apply_zero_boundary_conditions(spaces[0]->dof_map(), rhs_2);

		DVectord sol_c_1 = local_zeros(local_size(external_force));
		DVectord sol_c_2 = local_zeros(local_size(external_force));

		// DVectord rhs_c_1 = transpose(orthogonal_trafo) * transpose(transfer_operator) * rhs_1;
		// DVectord rhs_c_2 = transpose(orthogonal_trafo) * transpose(transfer_operator) * rhs_2;

		// DSMatrixd K_c_11  = transpose(orthogonal_trafo) *
		// DSMatrixd(
		// 	transpose(transfer_operator) * 
		// 	K_11 * 
		// 	transfer_operator) * 
		// orthogonal_trafo;


		// DSMatrixd K_c_21  = transpose(orthogonal_trafo) *
		// DSMatrixd(
		// 	transpose(transfer_operator) * 
		// 	K_21 * 
		// 	transfer_operator) * 
		// orthogonal_trafo;


		// DSMatrixd K_c_22  = transpose(orthogonal_trafo) *
		// DSMatrixd(
		// 	transpose(transfer_operator) * 
		// 	K_22 * 
		// 	transfer_operator) * 
		// orthogonal_trafo;


		//end: not contact
		DVectord &rhs_c_1 = rhs_1;
		DVectord &rhs_c_2 = rhs_2;

		DSMatrixd &K_11_c = K_11;
		DSMatrixd &K_21_c = K_21;
		DSMatrixd &K_22_c = K_22;
		//end: not contact

		DVectord rhs_c;
		DSMatrixd K_c;
		make_block_system(K_11_c, K_21_c, K_22_c, K_c, rhs_c_1, rhs_c_2, rhs_c);

		// SemismoothNewton<DSMatrixd, DVectord, PETSC_EXPERIMENTAL> newton(linear_solver);
		SemismoothNewton<DSMatrixd, DVectord> newton(linear_solver);
		newton.verbose(true);
		newton.max_it(40);

		DVectord dummy = local_values(local_size(rhs_1).get(0) + local_size(rhs_2).get(0), 1000000);
		newton.set_box_constraints(make_upper_bound_constraints(make_ref(dummy)));
		// newton.set_box_constraints(make_upper_bound_constraints(make_ref(gap)));

		DVectord sol_c = local_zeros(local_size(rhs_c));
		newton.solve(K_c, rhs_c, sol_c);


		unmake_block_system(sol_c, sol_c_1, sol_c_2);

		// displacement_increment = transfer_operator * (orthogonal_trafo * sol_c_1);		
		// total_displacement += displacement_increment;
		// velocity = orthogonal_trafo * sol_c_1;		

		//begin: not contact
		displacement_increment = sol_c_1;		
		total_displacement += displacement_increment;
		velocity = sol_c_1;	
		//end: not contact


		new_internal_force = stiffness_matrix * total_displacement;
		apply_zero_boundary_conditions(spaces[0]->dof_map(), new_internal_force);
	}

	void ContactProblem::contact_stabilized_newmark(const double dt)
	{
		//predictor step
		DVectord pred = orthogonal_trafo * utopia::min(orthogonal_trafo * (dt * velocity), gap);
		DVectord rhs  = internal_mass_matrix * pred + (dt*dt/2.) * (2. * external_force - internal_force);
		DSMatrixd K   = internal_mass_matrix + (dt*dt/4.) * stiffness_matrix;

		DVectord sol_c = local_zeros(local_size(external_force));
		DVectord rhs_c = transpose(orthogonal_trafo) * transpose(transfer_operator) * rhs;
		DSMatrixd K_c  = transpose(orthogonal_trafo) *
		DSMatrixd(
			transpose(transfer_operator) * 
			K * 
			transfer_operator) * 
		orthogonal_trafo;

		SemismoothNewton<DSMatrixd, DVectord> newton(linear_solver);
		newton.verbose(true);
		newton.max_it(40);

		newton.set_box_constraints(make_upper_bound_constraints(make_ref(gap)));
		newton.solve(K_c, rhs_c, sol_c);

		displacement_increment = transfer_operator * (orthogonal_trafo * sol_c);		
		total_displacement += displacement_increment;

		new_internal_force = stiffness_matrix * total_displacement;
		apply_zero_boundary_conditions(spaces[0]->dof_map(), new_internal_force);

		//compute acceleration
		DVectord dt_x_M_x_acc = (dt/2.) * (2. * external_force - internal_force - new_internal_force);
		DVectord vel_inc = local_zeros(local_size(dt_x_M_x_acc));
		apply_zero_boundary_conditions(spaces[0]->dof_map(), dt_x_M_x_acc);
		// solve(mass_matrix, dt_x_M_x_acc, vel_inc);
		vel_inc = e_mul(inverse_mass_vector, dt_x_M_x_acc);

		velocity += vel_inc;

		// velocity = 1/dt * displacement_increment;
	}

	void ContactProblem::step(const double dt)
	{

		if(verbose) std::cout << "Contact problem: computing step" << std::endl;
		
		disp("n_dofs:");
		disp(size(external_force));
		const int dim = mesh->mesh_dimension();

		if(iteration > 0) {
			Read<DVectord> r_od(old_displacement_increment), r_v(velocity), r_td(total_displacement);

			search_radius = 1e-2;
			auto r = range(old_displacement_increment);
			for(auto i = r.begin(); i < r.end(); i += dim) {
				double length = 0.0;
				double length_total_disp = 0.0;

				for(int j = 0; j < dim; ++j) {
					double v = old_displacement_increment.get(i + j) +  2. * dt * velocity.get(i + j);
					length += v * v;
					length_total_disp += total_displacement.get(i + j) *  total_displacement.get(i + j);
				}

				search_radius = std::max(search_radius, std::min(1.1 * std::sqrt(length), std::sqrt(length_total_disp)));
			}

			if(is_inpulse_) {
				external_force = zeros(size(external_force));
			}
		}

		search_radius += 1e-4;
		comm.all_reduce(&search_radius, 1, moonolith::MPIMax());
		disp("predicted search radius");
		disp(search_radius);

		compute_contact_conditions();

		MOONOLITH_EVENT_BEGIN("solving_vi");

		if(dynamic_contact) {
			// classic_newmark_with_contact(dt);
			// classic_newmark_with_contact_2(dt);
			contact_stabilized_newmark(dt);
			// contact_stabilized_newmark_monolithic(dt);
		} else {
			implicity_euler(dt);
		}
		
		old_displacement_increment = displacement_increment;
		internal_force = new_internal_force;

		compute_normal_stress(dt);
		apply_displacement(displacement_increment);
		compute_energy(dt);
		MOONOLITH_EVENT_END("solving_vi");
		++iteration;
	}


	void ContactProblem::compute_normal_stress(const double dt)
	{
		const int dim = mesh->mesh_dimension();
		normal_stress = local_zeros(local_size(total_displacement));
		DVectord unscaled_stress = boundary_mass_inv * (external_force - internal_force);
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
	}

	void ContactProblem::apply_displacement(const DVectord &displacement_increment)
	{	
		int sys_num = context_ptr->system.number();

		if(!comm.is_alone()) {
			auto r = range(displacement_increment);
			Read<DVectord> r_d(displacement_increment);

			auto m_begin = mesh->active_local_elements_begin();
			auto m_end   = mesh->active_local_elements_end();

			std::vector<PetscInt> idx;
			std::set<PetscInt> unique_idx;
			std::map<libMesh::dof_id_type, double> idx_to_value;
			std::vector<libMesh::dof_id_type> dof_indices;
			
			auto &dof_map = spaces[0]->dof_map();
			for(auto m_it = m_begin; m_it != m_end; ++m_it) { 
				dof_map.dof_indices(*m_it, dof_indices);
				for(auto dof_id : dof_indices) {
					if(r.inside(dof_id)) {
						idx_to_value[dof_id] = displacement_increment.get(dof_id);
					} else {
						unique_idx.insert(dof_id);
					}
				}
			}

			idx.insert(idx.end(), unique_idx.begin(), unique_idx.end());
			DVectord out = displacement_increment.select(idx);
			{
				Read<DVectord> r_out(out);
				auto range_out = range(out);

				for(std::size_t i = 0; i < idx.size(); ++i) {
					idx_to_value[idx[i]] = out.get(range_out.begin() + i);
				}
			}	

			for(auto m_it = m_begin; m_it != m_end; ++m_it) { 
				auto &e = **m_it;
				for(int i = 0; i < e.n_nodes(); ++i) {
					auto &node = e.node_ref(i);

					for(unsigned int c = 0; c < mesh->mesh_dimension(); ++c) {
						const int dof_id = node.dof_number(sys_num, c, 0);
						assert(idx_to_value.find(dof_id) != idx_to_value.end());
						double &val = idx_to_value[dof_id];
						node(c) += val;
						val = 0.;
					}
				}
			}

		} else {

			Read<DVectord> r_d(displacement_increment);

			auto m_it  = mesh->local_nodes_begin();
			auto m_end = mesh->local_nodes_end();


			for(; m_it != m_end; ++m_it) { 
				for(unsigned int c = 0; c < mesh->mesh_dimension(); ++c) {
					const int dof_id = (*m_it)->dof_number(sys_num, c, 0);
					(**m_it)(c) += displacement_increment.get(dof_id);
				}
			}
		}
	}

	void ContactProblem::save(const double dt, const std::string &output_dir)
	{

		const int dim = mesh->mesh_dimension();

		if(iteration > 0 && comm.is_alone() && utopia::Utopia::Instance().get("plot") == "true") {
			std::vector<double> is_contact_node_x(local_size(is_contact_node).get(0)/dim, 0.);

			{
				Read<DVectord> r_icn(is_contact_node);
				Range r = range(is_contact_node);
				for(std::size_t i = 0; i < is_contact_node_x.size(); ++i) {
					is_contact_node_x[i] = is_contact_node.get(r.begin() + i * dim);
				}
			}

			plot_mesh_f(*mesh, &is_contact_node_x[0], "time_series_m/m" + std::to_string(iteration));
		} 

		auto &aux = context_ptr->equation_systems.get_system<libMesh::ExplicitSystem>("aux");

		MeshBase::const_node_iterator nd = mesh->local_nodes_begin();
		const MeshBase::const_node_iterator nd_end = mesh->local_nodes_end();


		DVectord contact_stress = external_force - internal_force;
		DVectord u_contact_stress = local_zeros(local_size(contact_stress));
		solve(mass_matrix, contact_stress, u_contact_stress);
		apply_zero_boundary_conditions(spaces[0]->dof_map(), u_contact_stress);

		scale_normal_vector_with_gap(dim, normals, gap, rays);

		{
			Read<DVectord> r_d(total_displacement), r_v(velocity), r_f(internal_force), r_c(contact_stress);

			for (; nd != nd_end; ++nd)
			{
				const Node * node = *nd;
				for (unsigned int d = 0; d < dim; ++d) {
					unsigned int source_dof 	= node->dof_number(context_ptr->system.number(), d, 0);
					unsigned int dest_dof_disp 	= node->dof_number(aux.number(), var_num_aux[d], 0);
					unsigned int dest_dof_vel 	= node->dof_number(aux.number(), var_num_aux[dim+d], 0);
					unsigned int dest_dof_force = node->dof_number(aux.number(), var_num_aux[2*dim+d], 0);
					unsigned int dest_dof_contact_stress = node->dof_number(aux.number(), var_num_aux[3*dim+d], 0);

					aux.solution->set(dest_dof_disp, total_displacement.get(source_dof));
					aux.solution->set(dest_dof_vel, velocity.get(source_dof));
					aux.solution->set(dest_dof_force, internal_force.get(source_dof));
					aux.solution->set(dest_dof_contact_stress, u_contact_stress.get(source_dof));
				}

				unsigned int source_dof = node->dof_number(context_ptr->system.number(), 0, 0);
				unsigned int dest_dof_is_contact = node->dof_number(aux.number(), var_num_aux[4*dim], 0);
				aux.solution->set(dest_dof_is_contact, is_contact_node.get(source_dof));


				for (unsigned int d = 0; d < dim; ++d) {
					unsigned int source_dof = node->dof_number(context_ptr->system.number(), d, 0);
					unsigned int dest_dof_ray = node->dof_number(aux.number(), var_num_aux[4*dim + 1 + d], 0);
					aux.solution->set(dest_dof_ray, rays.get(source_dof));
				}
			}
		}

		convert(total_displacement, *context_ptr->system.solution);
		aux.solution->close();
		output->write_timestep(output_dir + "/sol_" + std::to_string(comm.size()) + ".e", context_ptr->equation_systems, iteration + 1, dt*(iteration + 1));
	}

	void ContactProblem::compute_energy(const double dt)
	{
		Energy e;
		e.t = dt * (iteration + 1);
		e.kinetic_energy = 0.5 * dot(mass_matrix * velocity, velocity);
		e.elastic_energy = 0.5 * dot(stiffness_matrix * total_displacement, total_displacement);
		e.potential_energy = e.elastic_energy - dot(external_force,  total_displacement);
		e.contact_energy = sum(mass_matrix * normal_stress);

		energy.push_back(e);
	}

	void ContactProblem::save_energy(const std::string &path)
	{
		if(!comm.is_root()) return;

		std::ofstream os(path.c_str());

		if(!os.good()) {
			std::cerr << "[Error] ContactProblem::save_energy: unable to write in file" << std::endl;
			return;
		}

		for(auto &e : energy) {
			os << e.t << " " << e.kinetic_energy << " "  << e.elastic_energy << " " << e.potential_energy << " " << e.contact_energy << "\n";
		}

		os.close();
	}
}
