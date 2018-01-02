#include "utopia_FSITest.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"
#include "utopia_assemble_volume_transfer.hpp"
#include "moonolith_communicator.hpp"

namespace utopia {

	static void apply_displacement(
		moonolith::Communicator &comm,
		const libMesh::DofMap &dof_map,
		const std::vector<int> &vars,
		const DVectord &disp,
		libMesh::MeshBase &mesh
		)
	{	
		int sys_num = 0;

		if(!comm.is_alone()) {
			auto r = range(disp);
			Read<DVectord> r_d(disp);

			auto m_begin = mesh.active_local_elements_begin();
			auto m_end   = mesh.active_local_elements_end();

			std::vector<PetscInt> idx;
			std::set<PetscInt> unique_idx;
			std::map<libMesh::dof_id_type, double> idx_to_value;
			std::vector<libMesh::dof_id_type> dof_indices;
			
			for(auto m_it = m_begin; m_it != m_end; ++m_it) { 
				dof_map.dof_indices(*m_it, dof_indices);
				for(auto dof_id : dof_indices) {
					if(r.inside(dof_id)) {
						idx_to_value[dof_id] = disp.get(dof_id);
					} else {
						unique_idx.insert(dof_id);
					}
				}
			}

			idx.insert(idx.end(), unique_idx.begin(), unique_idx.end());
			DVectord out = disp.select(idx);
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

					for(unsigned int c = 0; c < vars.size(); ++c) {
						const int dof_id = node.dof_number(sys_num, vars[c], 0);
						assert(idx_to_value.find(dof_id) != idx_to_value.end());
						double &val = idx_to_value[dof_id];
						node(c) += val;
						val = 0.;
					}
				}
			}

		} else {

			Read<DVectord> r_d(disp);

			auto m_it  = mesh.local_nodes_begin();
			auto m_end = mesh.local_nodes_end();


			for(; m_it != m_end; ++m_it) { 
				for(unsigned int c = 0; c < vars.size(); ++c) {
					const int dof_id = (*m_it)->dof_number(sys_num, vars[c], 0);
					(**m_it)(c) += disp.get(dof_id);
				}
			}
		}
	}

	void run_fsi_test(libMesh::LibMeshInit &init)
	{
		moonolith::Communicator comm(init.comm().get());

		// const unsigned int nx_fluid = 210;
		// const unsigned int ny_fluid = 70;
		// const unsigned int nx_solid = 20;
		// const unsigned int ny_solid = 80;

		// const unsigned int nx_fluid = 105;
		// const unsigned int ny_fluid = 35;
		// const unsigned int nx_solid = 10;
		// const unsigned int ny_solid = 40;

		const unsigned int nx_fluid = 50;
		const unsigned int ny_fluid = 15;
		const unsigned int nx_solid = 5;
		const unsigned int ny_solid = 20;

		////////////////////////////////////////////////////////////////////////////////////
		//Fluid discretization
		////////////////////////////////////////////////////////////////////////////////////

		auto fluid_mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());		
		libMesh::MeshTools::Generation::build_square(*fluid_mesh,
			nx_fluid, ny_fluid,
			-1.5, 1.5,
			0, 1.,
			libMesh::QUAD8);

		auto fluid_system = std::make_shared<libMesh::EquationSystems>(*fluid_mesh);	
		fluid_system->add_system<libMesh::LinearImplicitSystem>("fluid");

		auto V_fx = LibMeshFunctionSpace(fluid_system, libMesh::LAGRANGE, libMesh::SECOND, "vel_x");
		auto V_fy = LibMeshFunctionSpace(fluid_system, libMesh::LAGRANGE, libMesh::SECOND, "vel_y");
		auto V_f  = V_fx * V_fy;
		auto Q_f =  LibMeshFunctionSpace(fluid_system, libMesh::LAGRANGE, libMesh::FIRST, "pressure");

		auto u_f = trial(V_f);
		auto v_f = test(V_f);

		auto u_fx = u_f[0];
		auto u_fy = u_f[1];

		auto p_f = trial(Q_f);
		auto q_f = test(Q_f);

		////////////////////////////////////////////////////////////////////////////////////
		//Solid discretization
		////////////////////////////////////////////////////////////////////////////////////

		auto solid_mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());
		libMesh::MeshTools::Generation::build_square(*solid_mesh,
			nx_solid, ny_solid,
			-0.05, 0.05,
			0, .4,
			libMesh::QUAD4);

		auto solid_system = std::make_shared<libMesh::EquationSystems>(*solid_mesh);	
		solid_system->add_system<libMesh::LinearImplicitSystem>("solid");	

		
		auto V_sx = LibMeshFunctionSpace(solid_system, libMesh::LAGRANGE, libMesh::FIRST, "disp_x");
		auto V_sy = LibMeshFunctionSpace(solid_system, libMesh::LAGRANGE, libMesh::FIRST, "disp_y");
		auto V_s  = V_sx * V_sy;

		auto u_s = trial(V_s);
		auto v_s = test(V_s);

		auto u_sx = u_s[0];
		auto u_sy = u_s[1];

		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////

		const double mu_f  = 1.;
		const double rho_f = 0.1;
		const double mu_s  = 50.;
		const double lambda_s = 50.;

		const double dt = 0.01;
		
		const std::size_t n_ts = 100;

		DVectord sol_f;
		DVectord sol_fold;
		DVectord fsi_forcing_term_f;

		DVectord old_displacement_s;
		DVectord displacement_s;
		DVectord reaction_force_s;
		DVectord fsi_velocity_s;

		auto uk_f 	 = interpolate(sol_f, u_f);
		auto uk_fold = interpolate(sol_fold, u_f);
		auto fk_f    = interpolate(fsi_forcing_term_f, u_f);

		auto g_uk_f  = grad(uk_f);

		auto uk_s = interpolate(displacement_s, u_s);

		//equations fluid
		LMDenseVector one = values(v_f.codim(), 1.);
		auto mass_f     = integral(inner(u_f, v_f)) + integral(inner(p_f, q_f));
		auto mass_rhs_f = integral(inner(coeff(one), v_f)) + integral(inner(coeff(1.), q_f));

		auto b_form_f = integral(inner(u_f, v_f)) 
						+ dt * integral(inner(mu_f * (transpose(grad(u_f)) + grad(u_f)), grad(v_f))
						- (rho_f * dt) * inner(g_uk_f * u_f, v_f))
						+ dt * integral(inner(p_f, div(v_f)))
						+ integral(inner(div(u_f), q_f));

		auto l_form_f = integral(inner(coeff(0.), q_f))
						+ integral(inner(uk_fold, v_f))
						+ integral(dt * inner(fk_f, v_f)); 


		auto eq_fluid = equations( b_form_f == l_form_f );
		
		//equations solid
		auto e_u = 0.5 * ( transpose(grad(u_s)) + grad(u_s) ); 
		auto e_v = 0.5 * ( transpose(grad(v_s)) + grad(v_s) );
		LMDenseVector z = zeros(2);


		auto mass_s = integral(inner(u_s, v_s));
		auto mass_rhs_s = integral(inner(coeff(one), v_s));

		//force response
		auto eq_solid = integral((2. * mu_s) * inner(e_u, e_v) + lambda_s * inner(div(u_s), div(v_s))) == integral(inner(coeff(z), v_s));
		
		//constraints fluid
		auto constr_f = 
			constraints(
				boundary_conditions(u_fy == coeff(0.),   {0, 1, 2, 3}),
				boundary_conditions(u_fx == coeff(0.),   {0, 2}),
				boundary_conditions(u_fx == coeff(0.1),  {1, 3})
				);

		//constraints solid
	    auto constr_s =  constraints(
	    		boundary_conditions(u_sx == coeff(0.), {0}),
	    		boundary_conditions(u_sy == coeff(0.), {0})
	    	);


	    //FIXME should be hidden
	    init_constraints(constr_f);
	    init_constraints(constr_s);
	   
	    V_fx.initialize();
	    V_sx.initialize();
	    //

	    // sol_f    = local_zeros(V_fx.dof_map().n_local_dofs());
	    // sol_fold = local_zeros(V_fx.dof_map().n_local_dofs());

	    auto &dof_map = V_fx.dof_map();
	    dof_map.prepare_send_list();
	    sol_f    = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), dof_map.get_send_list());
	    sol_fold = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), dof_map.get_send_list());

	    fsi_forcing_term_f = local_zeros(V_fx.dof_map().n_local_dofs());
	    displacement_s = local_zeros(V_sx.dof_map().n_local_dofs());

	    DSMatrixd mat_s, mass_mat_s, mass_mat_f;
	    DVectord rhs_s,  mass_vec_s, mass_vec_f;

	    assemble_expression_v<LibMeshFunctionSpace>(eq_solid, mat_s, rhs_s);
	    assemble_expression_v<LibMeshFunctionSpace>(mass_f == mass_rhs_f, mass_mat_f,  mass_vec_f, false);
	    assemble_expression_v<LibMeshFunctionSpace>(mass_s == mass_rhs_s, mass_mat_s,  mass_vec_s, false);

	    NonLinearFEFunction<DSMatrixd, DVectord, decltype(eq_fluid)> nl_fun(eq_fluid);

	    auto linear_solver = std::make_shared<Factorization<DSMatrixd, DVectord>>();
	    Newton<DSMatrixd, DVectord> solver(linear_solver);
	    solver.verbose(true);
	    
	    libMesh::ExodusII_IO io_fluid(*fluid_mesh);
	    libMesh::ExodusII_IO io_solid(*solid_mesh);
	    	
	    for(std::size_t ts = 0; ts < n_ts; ++ts) {
	    	std::cout << "----------------------------\n";
	    	std::cout << "t: " << (ts * dt) << "/" << (n_ts * dt) << std::endl;
	    	//assemble transfer operator
	    	if(ts > 0) {
	    		apply_displacement(comm, V_sx.dof_map(), {0, 1}, displacement_s - old_displacement_s, *solid_mesh);
	    		// plot_mesh(*solid_mesh, "mesh/" + std::to_string(ts));
	    	}

	    	DSMatrixd B;
	    	
	    	if(!assemble_volume_transfer(
	    	    comm,
	    	    fluid_mesh,
	    	    solid_mesh,
	    	    make_ref(V_fx.dof_map()),
	    	    make_ref(V_sx.dof_map()),
	    	    V_fx.subspace_id(),
	    	    V_sx.subspace_id(),
	    	    false, 
	    	    V_s.n_subspaces(),
	    	    B))
	    	{

	    		std::cerr << "Volume transfer failure" << std::endl;
	    	}

 			fsi_velocity_s = local_zeros(local_size(displacement_s));
	    	old_displacement_s = displacement_s;
	    	sol_fold = sol_f;

	    	std::vector<libMesh::dof_id_type> pressure_index;
	    	Q_f.dof_map().local_variable_indices(pressure_index, Q_f.mesh(), Q_f.subspace_id());

	    	long pressure_dofs = pressure_index.size();
	    	comm.all_reduce(&pressure_dofs, 1, moonolith::MPISum());

	    	// Fixed point iteration for current time-step
	    	bool converged = false;
	    	int outer_iter = 0;
	    	while(!converged) {
	    		std::cout << "outer iter: " << outer_iter << std::endl;
	    		//transfer velocity to solid
	    		DVectord temp_s = B * sol_f;
	    		std::cout << "mag temp_s: " << double(norm2(temp_s)) << std::endl;
	    		linear_solver->solve(mass_mat_s, temp_s, fsi_velocity_s);
	    		std::cout << "mag vel: " << double(norm2(fsi_velocity_s)) << std::endl;

	    		DVectord displacement_prev = displacement_s;
	    		displacement_s = old_displacement_s + dt * fsi_velocity_s;

	    		const double diff = norm2(displacement_prev - displacement_s);
	    		std::cout << "diff: " << diff << std::endl;

	    		//FIXMe
	    		if(outer_iter > 0) {
	    			converged = diff < 1e-5;
	    			if(converged) break;
	    		}

	    		reaction_force_s = rhs_s - mat_s * displacement_s;
	    		std::cout << "mag_f: " << double(norm2(reaction_force_s)) << std::endl;

	    		//update the forcing term in the fluid
	    		linear_solver->solve(mass_mat_s, reaction_force_s, temp_s);
	    		DVectord temp_f = transpose(B) * temp_s;
	    		linear_solver->solve(mass_mat_f, temp_f, fsi_forcing_term_f);

	    		const double mag_fsi = norm2(fsi_forcing_term_f);
	    		std::cout << "mag_fsi: " << mag_fsi << std::endl;
	    		
	    		DVectord temp = sol_f;

	     		if(!solver.solve(nl_fun, sol_f)) {
	     			std::cerr << "FAILED TO SOLVE NONLINEAR SYSTEM" << std::endl;
	     			break;
	     		}

				comm.barrier();
	     		std::cout << "Here" << std::endl;

	     		double mean_pressure = 0.;
	     		
	     		{
	     			Read<DVectord> r(sol_f);
	     			for(auto idx : pressure_index) {
	     				mean_pressure += sol_f.get(idx);
	     			}	
	     		}

	     		comm.all_reduce(&mean_pressure, 1, moonolith::MPISum());
	     		mean_pressure /= pressure_dofs;

	     		DVectord sol_f_copy = sol_f;

	     		{
	     			Read<DVectord>  r(sol_f_copy);
	     			Write<DVectord> w(sol_f);
	     			for(auto idx : pressure_index) {
	     			 	sol_f.set(idx, sol_f_copy.get(idx) - mean_pressure);
	     			}	
	     		}
	     	
	     		// converged = true;
	     		++outer_iter;
	     		if(outer_iter > 10) break;
	     	}

	     	
	     	{
	     		//export fluid solution
	     		// convert(fsi_forcing_term_f, *V_fx.equation_system().solution);
	     		convert(sol_f, *V_fx.equation_system().solution);
	     		V_fx.equation_system().solution->close();
	     		io_fluid.write_timestep("fsi_fluid.e", V_fx.equation_systems(), ts + 1, ts * dt);
	     	}

	     	{
	     		//export solid solution
	     		convert(displacement_s, *V_sx.equation_system().solution);
	     		// convert(reaction_force_s, *V_sx.equation_system().solution);
	     		V_sx.equation_system().solution->close();
	     		io_solid.write_timestep("fsi_solid.e", V_sx.equation_systems(), ts + 1, ts * dt);
	     	}

	    }
	}
}