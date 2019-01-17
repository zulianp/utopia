#include "utopia_FractureFlowTransportSimulation.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_FractureFlowUtils.hpp"

#include <iostream>

namespace utopia {

	FractureFlowTransportSimulation::FractureFlowTransportSimulation(libMesh::Parallel::Communicator &comm)
	: steady_flow_(comm), preset_velocity_field_(false)
	{}

	void FractureFlowTransportSimulation::read(utopia::Input &in)
	{
		steady_flow_.read(in);

		auto &V_s = steady_flow_.matrix->space.subspace(0);
		auto &V_f = steady_flow_.fracture_newtork->space.subspace(0);

		transport_m_.set_steady_state_function_space(steady_flow_.matrix->space);
		transport_f_.set_steady_state_function_space(steady_flow_.fracture_newtork->space);

		in.get("transport", transport_m_);
		in.get("transport", transport_f_);

		//specialized input for matrix
		in.get("transport", [this](Input &in) {
			in.get("matrix", transport_m_);
		});

		in.get("preset-velocity-field", preset_velocity_field_);
	}

	void FractureFlowTransportSimulation::compute_transport()
	{
		transport_f_.init(steady_flow_.x_f, *steady_flow_.fracture_newtork);
		transport_m_.init(steady_flow_.x_m, *steady_flow_.matrix);

		if(preset_velocity_field_) {

			const SizeType dim = transport_m_.space->subspace(0).mesh().spatial_dimension();

			each_write(transport_m_.velocity, [dim](const SizeType i) -> double {
				const SizeType local_idx = i % (dim + 1);

				if(local_idx == 0) {
					return -10;
				}

				// if(local_idx < dim) {
				// 	return 0.;
				// }
				return 0.;
			});
		}

		//compute upwind matrix
		transport_f_.assemble_system(*steady_flow_.fracture_newtork);
		transport_m_.assemble_system(*steady_flow_.matrix);
		transport_m_.compute_upwind_operator();

		auto &V_m = transport_m_.space->space().last_subspace();

		UVector xkp1 = transport_m_.velocity;
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

		int n_timesteps = 1;
		for(double t = transport_m_.dt; t < transport_m_.simulation_time; t += transport_m_.dt) {
			
			transport_m_.add_mass(xkp1, rhs);
			rhs += transport_m_.dt * transport_m_.f;

			transport_m_.constrain_concentration(rhs);

			solver_m.apply(rhs, xkp1);

			utopia::convert(xkp1, *V_m.equation_system().solution);
			V_m.equation_system().solution->close();

			io_m.write_timestep("transient.e", V_m.equation_systems(), ++n_timesteps, t);
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

		c.stop();
		std::cout << "Overall time: " << c << std::endl;
	}

	void FractureFlowTransportSimulation::Transport::init(const UVector &pressure, FractureFlow &flow)
	{
		auto &V = flow.space.subspace(0);
		const int dim = V.mesh().spatial_dimension();

		assert(space);

		if(!space->initialized()) {

			space->set_space( utopia::make_unique<ProductFunctionSpace<LibMeshFunctionSpace>>() );
			auto &transport_system = V.equation_systems().add_system<libMesh::LinearImplicitSystem>("velocity");

			for(int i = 0; i < dim; ++i) {
				space->space() *= LibMeshFunctionSpace(transport_system, transport_system.add_variable("vel_" + std::to_string(i), libMesh::Order(V.order(0)), libMesh::LAGRANGE) );
			}

			space->space() *= LibMeshFunctionSpace(transport_system, transport_system.add_variable("c", libMesh::Order(V.order(0)), libMesh::LAGRANGE) );
			space->subspace(0).initialize();
		} 

		{
			auto &aux = V.equation_systems().add_system<libMesh::LinearImplicitSystem>("aux_2");
			aux_space *= LibMeshFunctionSpace(aux, aux.add_variable("R_a", libMesh::Order(V.order(0)), libMesh::LAGRANGE) );
			aux_space *= LibMeshFunctionSpace(aux, aux.add_variable("upwind", libMesh::Order(V.order(0)), libMesh::LAGRANGE) );
			aux_space.subspace(0).initialize();
		}

		assert(dim + 1 == int(space->space().n_subspaces()));
	
		auto W  = space->space().subspace(0, dim);
		auto &P = space->space().subspace(dim);

		auto p = trial(P);

		auto u = trial(W);
		auto v = test(W);


		copy_values(V, pressure, P, pressure_w);

		auto ph = interpolate(pressure_w, p);

		//FIXME If do not put paranthesis it gives priority to the minus instead of multiplication
		//and there is a bug with the unary -dot(a, b) apparently
		auto l_form = -(
			 	inner(
			 		flow.diffusion_tensor * grad(ph),
			 		ctx_fun(flow.sampler) * v
			 	) * dX
		);

		auto b_form = inner(trial(space->space()), test(space->space())) * dX;

		UVector M_x_v;
		utopia::assemble(l_form, M_x_v);
		utopia::assemble(b_form, mass_matrix);

		mass_matrix_inverse.update(make_ref(mass_matrix));

		velocity = local_zeros(local_size(M_x_v));

		mass_vector = sum(mass_matrix, 1);

		if(lump_mass_matrix) {	
			velocity = e_mul(M_x_v, 1./mass_vector);
		} else {
			mass_matrix_inverse.apply(M_x_v, velocity);
		}

		copy_values(P, pressure_w, P, velocity);
	}

	

	void FractureFlowTransportSimulation::Transport::update_output()
	{
		auto &V = space->subspace(0);
		auto &sys = V.equation_system();
		utopia::convert(velocity, *sys.solution);
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

	FractureFlowTransportSimulation::Transport::Transport() :
	lump_mass_matrix(false),
	simulation_time(1.),
	h1_regularization(false),
	regularization_parameter(1e-4),
	dt(0.1),
	use_upwinding(true)
	{}

	void FractureFlowTransportSimulation::Transport::read(Input &in)
	{
		in.get("lump-mass-matrix", lump_mass_matrix);
		in.get("dt", dt);
		in.get("simulation-time", simulation_time);
		in.get("h1-regularization", h1_regularization);
		in.get("regularization-parameter", regularization_parameter);
		in.get("use-upwinding", use_upwinding);

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
	}

	void FractureFlowTransportSimulation::Transport::compute_upwind_operator()
	{
		const int dim = space->subspace(0).mesh().spatial_dimension();

		auto W  = space->space().subspace(0, dim);
		auto &C = space->space().subspace(dim);

		std::cout << C.dof_map().n_variables() << std::endl;

		auto c = trial(C);
		auto q = test(C);
		auto v = trial(W);

		auto uh = interpolate(velocity, v);
		auto l_form = inner(uh, grad(q)) * dX;

		UVector R_a;
		utopia::assemble(l_form, R_a);

		transform_values(C, R_a, C, upwind_vector, [](const double &val) -> double {
			return (val >= 0)? 1. : -1.;
		});

		UVector R_aux;
		copy_values(C, R_a, aux_space.subspace(0), R_aux);
		copy_values(C, upwind_vector, aux_space.subspace(1), R_aux);


		utopia::convert(R_aux, *aux_space.subspace(0).equation_system().solution);
		aux_space.subspace(0).equation_system().solution->close();
	}

	void FractureFlowTransportSimulation::Transport::assemble_system(FractureFlow &flow)
	{
		const int dim = space->subspace(0).mesh().spatial_dimension();

		auto W  = space->space().subspace(0, dim);
		auto &C = space->space().subspace(dim);
		auto &mesh = C.mesh();
		auto &dof_map = C.dof_map();

		std::cout << C.dof_map().n_variables() << std::endl;

		auto c = trial(C);
		auto q = test(C);
		auto v = trial(W);

		auto uh = interpolate(velocity, v);
		auto ph = interpolate(pressure_w, c);

		auto sampler_fun = ctx_fun(flow.sampler);

		static const int NONE = 0, FULL = 1, ARTIFICIAL_DIFFUSION = 2, SUPG = 3;
		const int upwind_type = SUPG;


		//FIXME This bugs because the var ranges are found from uh instead of c
		// auto b_form = inner(c * uh, grad(q)) * dX;	

		// auto vel = sampler_fun * flow.diffusion_tensor * grad(ph);
		auto vel = sampler_fun * grad(ph);
		auto b_form = inner(c * vel,  grad(q)) * dX;
		
		
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
						auto artificial_diffusion = inner(-0.01 * grad(c), grad(q)) * dX;
						auto eval_diff = eval(artificial_diffusion, ctx);
						eval_b_form += eval_diff;
						break;
					}

					case FULL:
					{	
						auto R_a    = inner(flow.diffusion_tensor * grad(ph), 		sampler_fun * grad(q)) * dX;
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
			// remove_mass(M_x_c0, c0);
			c0 = e_mul(M_x_c0, 1./mass_vector);
		} else {
			c0 = local_zeros(C.dof_map().n_local_dofs());
		}


		copy_values(C, c0, C, velocity);

		f = local_zeros(local_size(velocity));

		for(auto tag : in_out_flow) {
			auto l_form = surface_integral(inner(vel * c, normal() * q), tag);

		}

		if(forcing_function) {
			UVector ff;
			forcing_function->eval(velocity, ff);
			f += ff;
			double norm_f = norm2(f);
			std::cout << "norm_f " << norm_f << std::endl;
		}
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

	}

}
