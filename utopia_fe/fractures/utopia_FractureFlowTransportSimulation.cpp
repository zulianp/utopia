#include "utopia_FractureFlowTransportSimulation.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_FractureFlowUtils.hpp"

#include <iostream>

namespace utopia {

	FractureFlowTransportSimulation::FractureFlowTransportSimulation(libMesh::Parallel::Communicator &comm)
	: steady_flow_(comm), preset_velocity_field_(false), transport_m_("matrix"), transport_f_("fracture_network")
	{}

	void FractureFlowTransportSimulation::read(utopia::Input &in)
	{
		steady_flow_.read(in);

		auto &V_s = steady_flow_.matrix->space.subspace(0);
		auto &V_f = steady_flow_.fracture_network->space.subspace(0);

		transport_m_.set_steady_state_function_space(steady_flow_.matrix->space);
		transport_f_.set_steady_state_function_space(steady_flow_.fracture_network->space);

		in.get("transport", transport_m_);
		in.get("transport", transport_f_);

		//specialized input for matrix
		in.get("transport", [this](Input &in) {
			in.get("matrix", transport_m_);
			in.get("fracture-network", transport_f_);
		});

		in.get("preset-velocity-field", preset_velocity_field_);
		in.get("transient-solve-strategy", transient_solve_strategy);
	}


	void FractureFlowTransportSimulation::write_result_csv(const USparseMatrix &A)
	{
		auto &V_m = transport_m_.space->space().last_subspace();
		auto &V_f = transport_f_.space->space().last_subspace();


		std::vector<std::string> header(4 + 2);

		header[0] = "n_cells(Omega_0)";
		header[1] = "n_cells(Omega_1)";
		header[2] = "n_cells(Omega_2)";
		header[3] = "n_cells(Omega_3)";
		header[4] = "n_dofs";
		header[5] = "n_nz";

		std::vector<double> stats(4 + 2, 0.);


		stats[2] = V_f.mesh().n_active_elem();
		stats[3] = V_m.mesh().n_active_elem();

		stats[4] = size(A).get(0);
		stats[5] = utopia::nnz(A, 0.);


		std::cout << "nnz: " << stats[5] << std::endl;
		assert(stats[5] > 0);

	
		const std::size_t n_flows_m = transport_m_.total_in_out_flow.size();
		const std::size_t n_flows_f = transport_f_.total_in_out_flow.size();


		for(std::size_t i = 0; i < n_flows_f; ++i) {
			header.push_back("flow(Omega_f)_" + std::to_string(transport_f_.in_out_flow[i]));
			stats.push_back(transport_f_.total_in_out_flow[i]);

			std::cout << ("flow(Omega_f)_" + std::to_string(transport_f_.in_out_flow[i])) << " " << transport_f_.total_in_out_flow[i] << std::endl;
		}

		for(std::size_t i = 0; i < n_flows_m; ++i) {
			header.push_back("flow(Omega_m)_" + std::to_string(transport_m_.in_out_flow[i]));
			stats.push_back(transport_m_.total_in_out_flow[i]);

			std::cout << ("flow(Omega_m)_" + std::to_string(transport_m_.in_out_flow[i])) << " " << transport_m_.total_in_out_flow[i] << std::endl;
		}

		CSVWriter csv;
		csv.open_file("result.csv");
		csv.write_table_row(header);
		csv.write_table_row(stats);
		csv.close_file();

	}

	void FractureFlowTransportSimulation::compute_transport_separate()
	{

		transport_f_.init(steady_flow_.x_f, *steady_flow_.fracture_network);
		transport_m_.init(steady_flow_.x_m, *steady_flow_.matrix);

		transport_f_.assemble_system(*steady_flow_.fracture_network);
		transport_m_.assemble_system(*steady_flow_.matrix);

		auto &V_m = transport_m_.space->space().last_subspace();

		UVector &xkp1 = transport_m_.concentration;
		UVector rhs = local_zeros(local_size(xkp1));

		apply_boundary_conditions(V_m.dof_map(), transport_m_.system_matrix, xkp1);
		
		libMesh::Nemesis_IO io_m(V_m.mesh());

		utopia::convert(xkp1, *V_m.equation_system().solution);

		V_m.equation_system().solution->close();
		io_m.write_timestep("transient.e", V_m.equation_systems(), 1, 0);

		Factorization<USparseMatrix, UVector> solver_m, solver_f;
		solver_m.update(make_ref(transport_m_.system_matrix));
		// solver_f.update(make_ref(transport_f_.system_matrix));

		// write("A.m", transport_m_.system_matrix);

		assert(!empty(transport_m_.f));


		transport_m_.post_process_time_step(0., *steady_flow_.matrix);


		int n_timesteps = 1;
		for(double t = transport_m_.dt; t < transport_m_.simulation_time; t += transport_m_.dt) {
			
			transport_m_.add_mass(xkp1, rhs);
			rhs += transport_m_.dt * transport_m_.f;

			transport_m_.constrain_concentration(rhs);

			solver_m.apply(rhs, xkp1);


			transport_m_.post_process_time_step(t, *steady_flow_.matrix);


			utopia::convert(xkp1, *V_m.equation_system().solution);
			V_m.equation_system().solution->close();

			io_m.write_timestep("transient.e", V_m.equation_systems(), ++n_timesteps, t);
		}

		// write_result_csv*
	}

	void FractureFlowTransportSimulation::compute_transport_monolithic()
	{
		transport_f_.init(steady_flow_.x_f, *steady_flow_.fracture_network);
		transport_m_.init(steady_flow_.x_m, *steady_flow_.matrix);

		transport_f_.assemble_system(*steady_flow_.fracture_network);
		transport_m_.assemble_system(*steady_flow_.matrix);

		auto &V_m = transport_m_.space->space().last_subspace();
		auto &V_f = transport_f_.space->space().last_subspace();

		UVector &x_m = transport_m_.concentration;
		UVector &x_f = transport_f_.concentration;

		UVector rhs_m = local_zeros(local_size(x_m));
		UVector rhs_f = local_zeros(local_size(x_f));

		USparseMatrix &A_m = transport_m_.system_matrix;
		USparseMatrix &A_f = transport_f_.system_matrix;

		apply_boundary_conditions(V_m.dof_map(), A_m, x_m);
		apply_boundary_conditions(V_f.dof_map(), A_f, x_f);

		USparseMatrix A = Blocks<USparseMatrix>(3, 3,
		{
		    make_ref(A_m),			  nullptr, 					make_ref(steady_flow_.B_t),
		    nullptr, 				  make_ref(A_f), 			make_ref(steady_flow_.D_t),
		    make_ref(steady_flow_.B), make_ref(steady_flow_.D), nullptr
		});

		// write("A.m", A);

		//lagrange mult
		UVector z = local_zeros(local_size(steady_flow_.B).get(0));
		UVector lagr = local_zeros(local_size(z));
	
		libMesh::Nemesis_IO io_m(V_m.mesh()), io_f(V_f.mesh());
		utopia::convert(x_m, *V_m.equation_system().solution);
		
		V_m.equation_system().solution->close();
		io_m.write_timestep("transient.e", V_m.equation_systems(), 1, 0);

		utopia::convert(x_f, *V_f.equation_system().solution);
		V_f.equation_system().solution->close();
		io_f.write_timestep("transient_f.e", V_f.equation_systems(), 1, 0);

		UVector x = blocks(x_m, x_f, lagr);


		// KSPSolver<USparseMatrix, UVector> op;
		// op.set_initial_guess_non_zero(false);
		// op.ksp_type(KSPPREONLY);
		// op.pc_type(PCLU);
		// op.solver_package(MATSOLVERSUPERLU_DIST);

		Factorization<USparseMatrix, UVector> op;
		op.describe(std::cout);
		
		
		op.update(make_ref(A));
		

		int n_timesteps = 1;
		const double dt = transport_m_.dt;
		const double simulation_time = transport_m_.simulation_time;



		transport_f_.post_process_time_step(0., *steady_flow_.fracture_network);
		transport_m_.post_process_time_step(0., *steady_flow_.matrix);

		for(double t = dt; t < simulation_time; t += dt) {
			transport_m_.add_mass(x_m, rhs_m);
			transport_f_.add_mass(x_f, rhs_f);

			rhs_m += transport_m_.dt * transport_m_.f;
			rhs_f += transport_f_.dt * transport_f_.f;

			transport_m_.constrain_concentration(rhs_m);
			transport_f_.constrain_concentration(rhs_f);

			UVector rhs = blocks(rhs_m, rhs_f, z);
			
			op.apply(rhs, x);

			undo_blocks(x, x_m, x_f, lagr);


			transport_f_.post_process_time_step(t, *steady_flow_.fracture_network);
			transport_m_.post_process_time_step(t, *steady_flow_.matrix);

			utopia::convert(x_m, *V_m.equation_system().solution);
			V_m.equation_system().solution->close();

			utopia::convert(x_f, *V_f.equation_system().solution);
			V_f.equation_system().solution->close();

			io_m.write_timestep("transient.e",   V_m.equation_systems(), ++n_timesteps, t);
			io_f.write_timestep("transient_f.e", V_f.equation_systems(), ++n_timesteps, t);
		}

		write_result_csv(A);
	}


	void FractureFlowTransportSimulation::compute_transport()
	{

		if(steady_flow_.solve_strategy == "separate" || transient_solve_strategy == "separate") {
			compute_transport_separate();
		} else {
			compute_transport_monolithic();
		}
	}

	void FractureFlowTransportSimulation::write_output()
	{
		transport_m_.update_output();
		transport_f_.update_output();
		steady_flow_.write_output();
	}

	void FractureFlowTransportSimulation::run()
	{
		Chrono c;
		c.start();

		steady_flow_.assemble_systems();
		steady_flow_.init_coupling();

		if(steady_flow_.solve()) {
			compute_transport();
			write_output();
		}

		transport_f_.finalize();
		transport_m_.finalize();

		c.stop();
		std::cout << "Overall time: " << c << std::endl;
	}

	void FractureFlowTransportSimulation::Transport::assemble_aux_quantities(FractureFlow &flow)
	{
		//gradient recovery using the standard L2-projection

		m_utopia_status("FractureFlowTransportSimulation::Transport::assemble_aux_quantities begin");
		Chrono chrono; chrono.start();

		auto &C = space->subspace(0);
		const int dim = C.mesh().spatial_dimension();

		auto &P = aux_space.subspace(0);
		auto W = aux_space.subspace(1, dim + 1);

		auto p = trial(P);
		auto u = trial(W);
		auto v = test(W);

		UVector aux_pressure = ghosted(P.dof_map().n_local_dofs(), P.dof_map().n_dofs(), P.dof_map().get_send_list()); 
		copy_values(C, pressure_w, P, aux_pressure);
		synchronize(aux_pressure);

		auto ph = interpolate(aux_pressure, p);

		//FIXME If do not put paranthesis it gives priority to the minus instead of multiplication
		//and there is a bug with the unary -dot(a, b) apparently
		auto l_form = -(
			 	inner(
			 		flow.diffusion_tensor * grad(ph),
			 		ctx_fun(flow.sampler) * v
			 	) * dX
		);

		auto b_form = inner(trial(aux_space), test(aux_space)) * dX;

		USparseMatrix aux_mass_matrix;
		UVector M_x_v;
		utopia::assemble(l_form, M_x_v);
		utopia::assemble(b_form, aux_mass_matrix);

		UVector aux_values = local_zeros(local_size(M_x_v));

		//this mass vector contains also the porosity
		copy_values(C, mass_vector, aux_space.subspace(dim + 1), M_x_v);

		if(lump_mass_matrix) {	
			UVector aux_mass_vector = sum(aux_mass_matrix, 1);
			aux_values = e_mul(M_x_v, 1./aux_mass_vector);
		} else {
			GMRES<USparseMatrix, UVector>("bjacobi").solve(aux_mass_matrix, M_x_v, aux_values);
		}

		copy_values(C, pressure_w, P, aux_values);

		utopia::convert(aux_values, *P.equation_system().solution);
		P.equation_system().solution->close();


		m_utopia_status("FractureFlowTransportSimulation::Transport::assemble_aux_quantities end");
		chrono.stop();

		std::cout << chrono << std::endl;
	}

	void FractureFlowTransportSimulation::Transport::finalize()
	{
		csv.close_file();
	}

	void FractureFlowTransportSimulation::Transport::init(const UVector &pressure, FractureFlow &flow)
	{
		auto &V = flow.space.subspace(0);
		const int dim = V.mesh().spatial_dimension();

		assert(space);

		{
			csv.open_file(name + ".csv");
			std::vector<std::string> vals = { "t" };

			for(auto io : in_out_flow) {
				vals.push_back("in/out-flow" + std::to_string(io));
			}

			csv.write_table_row(vals);
		}

		if(!space->initialized()) {
			space->set_space( utopia::make_unique<ProductFunctionSpace<LibMeshFunctionSpace>>() );
			auto &transport_system = V.equation_systems().add_system<libMesh::LinearImplicitSystem>("concentration");
			space->space() *= LibMeshFunctionSpace(transport_system, transport_system.add_variable("c", libMesh::Order(V.order(0)), libMesh::LAGRANGE) );
			space->subspace(0).initialize();
		} 

		auto &C = space->space().subspace(0);

		pressure_w = ghosted(C.dof_map().n_local_dofs(), C.dof_map().n_dofs(), C.dof_map().get_send_list()); 
		concentration = ghosted(C.dof_map().n_local_dofs(), C.dof_map().n_dofs(), C.dof_map().get_send_list()); 

		copy_values(V, pressure, C, pressure_w);


		synchronize(pressure_w);

		{
			auto &aux = V.equation_systems().add_system<libMesh::LinearImplicitSystem>("aux_2");
			aux_space *= LibMeshFunctionSpace(aux, aux.add_variable("p", libMesh::Order(V.order(0)), libMesh::LAGRANGE) );
			
			for(int i = 0; i < dim; ++i) {
				aux_space *= LibMeshFunctionSpace(aux, aux.add_variable("vel_" + std::to_string(i), libMesh::Order(V.order(0)), libMesh::LAGRANGE) );
			}

			aux_space *= LibMeshFunctionSpace(aux, aux.add_variable("porosity", libMesh::Order(V.order(0)), libMesh::LAGRANGE) );

			aux_space.subspace(0).initialize();
		}

		assert(1 == int(space->space().n_subspaces()));
			
		auto b_form = inner(ctx_fun(porosity) * trial(space->space()), test(space->space())) * dX;
		// auto b_form = inner(trial(space->space()), test(space->space())) * dX;

		utopia::assemble(b_form, mass_matrix);

		mass_vector = sum(mass_matrix, 1);
		mass_matrix_inverse.update(make_ref(mass_matrix));

		double vol_x_porsity = sum(mass_vector);
		std::cout << "vol_x_porsity: " << vol_x_porsity << std::endl;

		assemble_aux_quantities(flow);
	}

	
	void FractureFlowTransportSimulation::Transport::update_output()
	{
		auto &V = space->subspace(0);
		auto &sys = V.equation_system();
		utopia::convert(concentration, *sys.solution);
		sys.solution->close();
	}

	void FractureFlowTransportSimulation::Transport::remove_mass(const UVector &in, UVector &out)
	{
		if(empty(out)) {
			out = local_zeros(local_size(in));
		}

		if(lump_mass_matrix) {
			out = e_mul(in, 1./mass_vector);
		} else {
			const bool ok = mass_matrix_inverse.apply(in, out); assert(ok); (void) ok;
		}
	}

	void FractureFlowTransportSimulation::Transport::add_mass(const UVector &in, UVector &out) const
	{
		if(empty(out)) {
			out = local_zeros(local_size(in));
		}

		if(lump_mass_matrix) {
			out = e_mul(in, mass_vector);
		} else {
			out = mass_matrix * in;
		}
	}

	FractureFlowTransportSimulation::Transport::Transport(const std::string &name) :
	name(name),
	lump_mass_matrix(false),
	simulation_time(1.),
	h1_regularization(false),
	regularization_parameter(0.5),
	dt(0.1),
	use_upwinding(true),
	boundary_factor(1.)
	{}

	void FractureFlowTransportSimulation::Transport::read(Input &in)
	{
		in.get("lump-mass-matrix", lump_mass_matrix);
		in.get("dt", dt);
		in.get("simulation-time", simulation_time);
		in.get("h1-regularization", h1_regularization);
		in.get("regularization-parameter", regularization_parameter);
		in.get("use-upwinding", use_upwinding);
		in.get("boundary-factor", boundary_factor);


		in.get("outflow", 
			[this](Input &in) {
				in.get_all([this](Input &in) {
					int side = -1;
					in.get("side", side);

					if(side >= 0) {
						in_out_flow.push_back(side);
					}
				});
		});

		in.get("inflow", 
			[this](Input &in) {
				in.get_all([this](Input &in) {
					int side = -1;
					in.get("side", side);

					if(side >= 0) {
						in_out_flow.push_back(side);
					}
				});
		});

		in.get("box", [this](Input &in) {
			box_min.resize(3, -std::numeric_limits<double>::max());
			box_max.resize(3,  std::numeric_limits<double>::max());

			in.get("x-min", box_min[0]);
			in.get("y-min", box_min[1]);
			in.get("z-min", box_min[2]);

			in.get("x-max", box_max[0]);
			in.get("y-max", box_max[1]);
			in.get("z-max", box_max[2]);
		});

		in.get("space", *space);

		if(space->initialized()) {
			forcing_function = utopia::make_unique< UIForcingFunction<LibMeshFunctionSpace, UVector> >(space->space().last_subspace());
			in.get("forcing-function", *forcing_function);
		} else {
			std::cerr << "[Warning] no space and forcing-function applied for simulation" << std::endl;
		}

		in.get("porosity", [this](Input &in) {
			auto subdomain_fun = utopia::make_unique<UISubdomainFunction<double>>();
			subdomain_fun->read(in);
			
			if(subdomain_fun->good()) {
				
				if(!subdomain_fun->has_default()) {
					subdomain_fun->set_default(utopia::make_unique<UIConstantFunction<double>>(1.));
				}

				porosity = std::move(subdomain_fun);
			} else {
				std::cerr << "[Error] improper porosity format" << std::endl;
			}



		});

		if(!porosity) {
			porosity = utopia::make_unique<UIConstantFunction<double>>(1.);
		}
	}

	void FractureFlowTransportSimulation::compute_upwind_operator()
	{


	}

	void FractureFlowTransportSimulation::Transport::assemble_system(FractureFlow &flow)
	{
		const int dim = space->subspace(0).mesh().spatial_dimension();

		auto &C = space->space().subspace(0);
		auto &mesh = C.mesh();
		auto &dof_map = C.dof_map();

		std::cout << C.dof_map().n_variables() << std::endl;

		auto c = trial(C);
		auto q = test(C);
		auto ph = interpolate(pressure_w, c);

		auto sampler_fun = ctx_fun(flow.sampler);

		static const int NONE = 0, FULL = 1, ARTIFICIAL_DIFFUSION = 2, SUPG = 3;
		const int upwind_type = SUPG;


		//FIXME This bugs because the var ranges are found from uh instead of c
		// auto b_form = inner(c * uh, grad(q)) * dX;	

		auto vel = sampler_fun * flow.diffusion_tensor * grad(ph);
		// auto b_form = inner(c * vel,  grad(q)) * dX;
		auto b_form = (inner(inner(-grad(c), vel), q) * dX);
		
		
		if(!use_upwinding) {
			utopia::assemble(b_form, gradient_matrix);
		} else {

			auto n_local_dofs = dof_map.n_local_dofs();
			gradient_matrix = local_sparse(n_local_dofs, n_local_dofs, max_nnz_x_row(C));

			Write<USparseMatrix> w_(gradient_matrix);

			std::vector<libMesh::dof_id_type> dofs;
			AssemblyContext<LIBMESH_TAG> ctx;
			for(auto e_it = elements_begin(mesh); e_it != elements_end(mesh); ++e_it) {
				ctx.set_current_element((*e_it)->id());
				ctx.set_has_assembled(false);
				ctx.init( b_form );

				dof_map.dof_indices(*e_it, dofs);

				auto eval_b_form = eval(b_form, ctx);

				//upwinding begin
				
				auto n = dofs.size();

				switch(upwind_type) {
					case SUPG:
					{	
						auto C = regularization_parameter/norm2(vel);

						auto supg = 
						(
							inner(
								inner(vel, grad(c)),
								inner(vel, grad(q)) * C
								) 
						) * dX;

						//auto hmin = (*e_it)->hmin();
						auto hmax = (*e_it)->hmax();
						auto eval_supg = eval(supg, ctx);
						
						eval_supg *= hmax;
						eval_b_form += eval_supg;
						break;
					}
					case ARTIFICIAL_DIFFUSION:
					{
						auto artificial_diffusion = inner(regularization_parameter * grad(c), grad(q)) * dX;
						auto eval_diff = eval(artificial_diffusion, ctx);
						eval_b_form += eval_diff;
						break;
					}

					case FULL:
					{	
						auto R_a       = inner(flow.diffusion_tensor * grad(ph), sampler_fun * grad(q)) * dX;
						auto eval_R_a  = eval(R_a, ctx);

						//DOES not work
						std::vector<bool> upwind_node(n, false);
						std::vector<double> d_total_out(n, 0.);

						double total_in = 0., total_out = 0.;

						for(std::size_t i = 0; i < n; ++i) {
							upwind_node[i] = eval_R_a.get(i) >= 0.;

							if(upwind_node[i]) {
								eval_b_form.set(i, i, eval_b_form.get(i, i) + eval_R_a.get(i));
								total_out += eval_R_a.get(i);
								d_total_out[i] += eval_R_a.get(i);
							} else {
								total_in -= eval_R_a.get(i);
							}
						}

						for(std::size_t i = 0; i < n; ++i) {
							if(!upwind_node[i]) {
								auto val = eval_b_form.get(i, i);
								eval_b_form.set(i, i, val * d_total_out[i]/total_in);
							}	

						}

						break;
					}

					default:
					{
						assert(false);
						break;
					}
				}	

				//upwinding end
				add_matrix(utopia::raw_type(eval_b_form), dofs, dofs, gradient_matrix);
			}

		}


		for(auto tag : in_out_flow) {
			auto flow_form = surface_integral(inner(vel * c, normal() * q), tag);

			auto temp_boundary_flow_matrix = std::make_shared<USparseMatrix>();
			partial_boundary_flow_matrix.push_back(temp_boundary_flow_matrix);

			utopia::assemble(flow_form, *temp_boundary_flow_matrix);
			if(empty(boundary_flow_matrix)) {
				boundary_flow_matrix = *temp_boundary_flow_matrix;
			} else {
				boundary_flow_matrix += *temp_boundary_flow_matrix;
			}

			std::cout << "boundary flow at " << tag << std::endl;
		}

		if(!empty(boundary_flow_matrix)) {
			gradient_matrix -= boundary_factor * boundary_flow_matrix;
		}

		if(lump_mass_matrix) {
			system_matrix = dt * gradient_matrix;
			system_matrix += USparseMatrix(diag(mass_vector));
		} else {
			system_matrix = mass_matrix + dt * gradient_matrix;
		}

		if(h1_regularization) {
			auto lapl = inner(grad(c), grad(q)) * dX;

			USparseMatrix lapl_matrix;
			utopia::assemble(lapl, lapl_matrix);

			system_matrix += (regularization_parameter * dt) * lapl_matrix;
		}

		UVector c0;

		if(!box_min.empty() && !box_max.empty()) {
			auto IV = boxed_fun(box_min, box_max,
				[](const std::vector<double> &x) -> double {
					return 10.;
				}
			);

			auto l_form = inner(ctx_fun(IV), q) * dX;

			UVector M_x_c0;
			utopia::assemble(l_form, M_x_c0);
			c0 = e_mul(M_x_c0, 1./mass_vector);
		} else {
			c0 = local_zeros(C.dof_map().n_local_dofs());
		}

		copy_values(C, c0, C, concentration);
		synchronize(concentration);

		f = local_zeros(local_size(concentration));
		
		UVector ff;

		if(forcing_function) {
			forcing_function->eval(concentration, ff);
			f += ff;
			double norm_f = norm2(f);
			std::cout << "norm_f " << norm_f << std::endl;
		}
	}

	void FractureFlowTransportSimulation::Transport::post_process_time_step(const double t, FractureFlow &flow)
	{
		std::vector<double> vals = { t };

		if(total_in_out_flow.empty()){
			total_in_out_flow.resize(partial_boundary_flow_matrix.size(), 0.);
		}

		std::size_t idx = 0;
		for(auto pbfm : partial_boundary_flow_matrix) {
			double outflow_val = sum(*pbfm * concentration);
			vals.push_back(outflow_val);

			total_in_out_flow[idx++] += outflow_val;

			std::cout << outflow_val << std::endl;
		}

		csv.write_table_row(vals);
	}

	void FractureFlowTransportSimulation::Transport::constrain_concentration(UVector &vec)
	{	
		auto &V = space->space().last_subspace();
		auto &dof_map = V.dof_map();

		const int dim = V.mesh().spatial_dimension();

		const bool has_constaints = dof_map.constraint_rows_begin() != dof_map.constraint_rows_end();

		const libMesh::DofConstraintValueMap &rhs_values = dof_map.get_primal_constraint_values();

		{
			Write<UVector> w_v(vec);

			if(has_constaints) {
				Range r = range(vec);
				for(SizeType i = r.begin(); i < r.end(); ++i) {
					
					if(dof_map.is_constrained_dof(i)) {

						auto valpos = rhs_values.find(i);

						// if(valpos != rhs_values.end()) {
						vec.set(i, (valpos == rhs_values.end()) ? 0 : valpos->second);
						// }
					}
				}
			}
		}

		synchronize(vec);

	}

}
