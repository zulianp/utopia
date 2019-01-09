#include "utopia_FractureFlowTransportSimulation.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_FractureFlowUtils.hpp"

#include <iostream>

namespace utopia {

	FractureFlowTransportSimulation::FractureFlowTransportSimulation(libMesh::Parallel::Communicator &comm)
	: steady_flow_(comm)
	{}

	void FractureFlowTransportSimulation::read(utopia::Input &in)
	{
		steady_flow_.read(in);

		auto &V_s = steady_flow_.matrix->space.subspace(0);
		auto &V_f = steady_flow_.fracture_newtork->space.subspace(0);

		transport_m_.lump_mass_matrix = false;
		transport_f_.lump_mass_matrix = false;

		in.get("lump-mass-matrix", transport_m_.lump_mass_matrix);
		in.get("lump-mass-matrix", transport_f_.lump_mass_matrix);


		simulation_time_ = 1.;
		in.get("simulation-time", simulation_time_);

		transport_f_.dt = 0.1;
		transport_m_.dt = 0.1;

		in.get("dt", transport_f_.dt);
		in.get("dt", transport_m_.dt);


	}

	void FractureFlowTransportSimulation::compute_transport()
	{
		transport_f_.init(steady_flow_.x_f, *steady_flow_.fracture_newtork);
		transport_m_.init(steady_flow_.x_m, *steady_flow_.matrix);

		//compute upwind matrix
		transport_f_.assemble_system();
		transport_f_.assemble_system();


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
		auto &aux = V.equation_systems().add_system<libMesh::LinearImplicitSystem>("velocity");

		space = utopia::make_unique<ProductFunctionSpace<LibMeshFunctionSpace>>();

		const int dim = V.mesh().spatial_dimension();

		for(int i = 0; i < dim; ++i) {
			(*space) *= LibMeshFunctionSpace(aux, aux.add_variable("vel_" + std::to_string(i), libMesh::Order(V.order(0)), libMesh::LAGRANGE) );
		}

		(*space) *= LibMeshFunctionSpace(aux, aux.add_variable("p", libMesh::Order(V.order(0)), libMesh::LAGRANGE) );

		space->subspace(0).initialize();

		auto W  = space->subspace(0, dim);
		auto &P = space->subspace(dim);

		auto p = trial(P);

		auto u = trial(W);
		auto v = test(W);

		UVector pressure_w;
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

		auto b_form = inner(trial(*space), test(*space)) * dX;

		UVector M_x_v;
		utopia::assemble(l_form, M_x_v);
		utopia::assemble(b_form, mass_matrix);

		velocity = local_zeros(local_size(M_x_v));

		if(lump_mass_matrix) {
			mass_vector = sum(mass_matrix, 1);
			velocity = e_mul(M_x_v, 1./mass_vector);
		} else {
			utopia::solve(mass_matrix, M_x_v, velocity);
		}

		copy_values(P, pressure_w, P, velocity);
	}

	void FractureFlowTransportSimulation::Transport::assemble_system()
	{
		const int dim = space->subspace(0).mesh().spatial_dimension();

		auto W  = space->subspace(0, dim);
		auto &C = space->subspace(dim);

		std::cout << C.dof_map().n_variables() << std::endl;

		auto c = trial(C);
		auto q = test(C);
		auto v = trial(W);

		auto uh = interpolate(velocity, v);

		//FIXME This bugs because the ranges are found from uh instead of c
		// auto b_form = inner(c * uh, grad(q)) * dX;	
		
		auto b_form = inner(c * uh, grad(q)) * dX;
		utopia::assemble(b_form, gradient_matrix);

		system_matrix = mass_matrix + dt * gradient_matrix;
	}

	void FractureFlowTransportSimulation::Transport::update_output()
	{
		auto &V = space->subspace(0);
		auto &sys = V.equation_system();
		utopia::convert(velocity, *sys.solution);
		sys.solution->close();
	}

}
