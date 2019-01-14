#include "utopia_SteadyFractureFlowSimulation.hpp"

#include "utopia.hpp"
#include "utopia_Socket.hpp"
#include "utopia_FractureFlowUtils.hpp"
#include "utopia_TransferUtils.hpp"
#include "utopia_TransferAssembler.hpp"

#include <iostream>

namespace utopia {

	SteadyFractureFlowSimulation::SteadyFractureFlowSimulation(libMesh::Parallel::Communicator &comm)
	{
		matrix 				= utopia::make_unique<FractureFlow>(comm);
		fracture_newtork 	= utopia::make_unique<FractureFlow>(comm);
		lagrange_multiplier = utopia::make_unique<FractureFlow>(comm);

		solve_strategy = "monolithic";
		use_mg = false;
		mg_sweeps = 1;
		mg_levels = 5;
		plot_matrix = false;
		write_operators_to_disk = false;
	}

	void SteadyFractureFlowSimulation::read(Input &is)
	{
		Chrono c;
		c.start();

		is.get("master", 	 *matrix);
		is.get("slave",  	 *fracture_newtork);
		is.get("multiplier", *lagrange_multiplier);

		is.get("solve-strategy", solve_strategy);
		is.get("use-mg", use_mg);
		is.get("mg-sweeps", mg_sweeps);
		is.get("mg-levels", mg_levels);
		is.get("plot-matrix", plot_matrix);

		if(plot_matrix) {
		    plot_mesh(matrix->mesh.mesh(), "matrix");
		}

		matrix->describe();
		fracture_newtork->describe();
		lagrange_multiplier->describe();

		std::cout << "solve-strategy: "  << solve_strategy << std::endl;
		std::cout << "use-mg:         "  << use_mg         << std::endl;
		if(use_mg) {
		    std::cout << "mg-sweeps:      "  << mg_sweeps      << std::endl;
		    std::cout << "mg-levels:      "  << mg_levels      << std::endl;
		}

		auto &V_m = matrix->space.subspace(0);
		auto &V_f = fracture_newtork->space.subspace(0);

		auto aux_matrix = utopia::make_unique<FractureFlowAuxSystem>(V_m);
		aux_matrix->sample(matrix->sampler);

		auto aux_fracture_newtork = utopia::make_unique<FractureFlowAuxSystem>(V_f);
		aux_fracture_newtork->sample(fracture_newtork->sampler);

		std::cout << "n_dofs: " << V_m.dof_map().n_dofs() << " + " <<  V_f.dof_map().n_dofs() << " + ";
		
		if(lagrange_multiplier->empty()) {
		    std::cout << V_f.dof_map().n_dofs();
		} else {
		    //Init mult space
		    lagrange_multiplier->space.subspace(0).initialize();
		    std::cout << lagrange_multiplier->space.subspace(0).dof_map().n_dofs();
		}

		std::cout << std::endl;

		c.stop();
		std::cout << "set-up time: " << c << std::endl;
	}

	void SteadyFractureFlowSimulation::assemble_systems()
	{
		Chrono c;
		c.start();

		auto &V_m = matrix->space.subspace(0);
		auto &V_f = fracture_newtork->space.subspace(0);

		auto u_m = trial(V_m);
		auto v_m = test(V_m);

		auto u_s = trial(V_f);
		auto v_f = test(V_f);

		auto eq_m = inner(matrix->diffusion_tensor * grad(u_m), ctx_fun(matrix->sampler) * grad(v_m)) * dX;
		auto eq_f = inner(fracture_newtork->diffusion_tensor  * grad(u_s), ctx_fun(fracture_newtork->sampler)  * grad(v_f)) * dX;

		utopia::assemble(eq_m, A_m);
		utopia::assemble(eq_f, A_f);

		x_m = local_zeros(V_m.dof_map().n_local_dofs());
		x_f = local_zeros(V_f.dof_map().n_local_dofs());

		matrix->forcing_function->eval(x_m, rhs_m);
		fracture_newtork->forcing_function->eval(x_f, rhs_f);

		apply_boundary_conditions(V_m.dof_map(), A_m, rhs_m);
		apply_boundary_conditions(V_f.dof_map(), A_f, rhs_f);


		c.stop();
		std::cout << "model assemly time: " << c << std::endl;

		if(write_operators_to_disk) {
		    write("A_m.m", A_m);
		    write("A_f.m", A_f);
		}
	}

	void SteadyFractureFlowSimulation::init_coupling()
	{
		Chrono c;
		c.start();

		auto &V_m = matrix->space.subspace(0);
		auto &V_f = fracture_newtork->space.subspace(0);

		if(lagrange_multiplier->empty()) {
			assemble_projection(V_m, V_f, B, D);
		} else {
			auto &V_l = lagrange_multiplier->space.subspace(0);
			assemble_projection(V_m, V_f, V_l, B, D);

		}

		D *= -1.;

		D_t = transpose(D);
		B_t = transpose(B);

		set_zero_at_constraint_rows(V_m.dof_map(), B_t);
		set_zero_at_constraint_rows(V_f.dof_map(), D_t);

		c.stop();
		std::cout << "transfer assemly time: " << c << std::endl;
	}

	bool SteadyFractureFlowSimulation::solve()
	{
		Chrono c;
		c.start();

		bool ok = false;

		if(solve_strategy == "staggered") {
 			ok = solve_cg_dual();
		} else if(solve_strategy == "separate") {
			ok = solve_separate();
		} else {
			ok = solve_monolithic();
		}

		c.stop();
		std::cout << "Solver time: " << c << std::endl;
		return ok;
	}

	bool SteadyFractureFlowSimulation::solve_cg_dual()
	{
		auto &V_m = matrix->space.subspace(0);
		auto &V_f = fracture_newtork->space.subspace(0);

		SPBlockConjugateGradient<USparseMatrix, UVector> solver;
		solver.verbose(true);
		solver.max_it(2000);
		solver.atol(1e-14);

		solver.use_simple_preconditioner();

		if(use_mg) {
		    auto mg = make_mg_solver(V_m, mg_levels);
		    solver.set_master_solver(mg);

		    solver.set_master_sweeps(mg_sweeps);
		    solver.set_master_max_it(mg->max_it());
		}

		solver.update(
		    make_ref(A_m),
		    make_ref(A_f),
		    make_ref(B),
		    make_ref(D),
		    make_ref(B_t),
		    make_ref(D_t)
		);

		return solver.apply(rhs_m, rhs_f, x_m, x_f, lagr);
	}

	bool SteadyFractureFlowSimulation::solve_monolithic()
	{
		USparseMatrix A = Blocks<USparseMatrix>(3, 3,
		{
		    make_ref(A_m), nullptr, make_ref(B_t),
		    nullptr, make_ref(A_f), make_ref(D_t),
		    make_ref(B), make_ref(D), nullptr
		});

		UVector z = local_zeros(local_size(rhs_f));
		UVector rhs = blocks(rhs_m, rhs_f, z);

		x_m = local_zeros(local_size(rhs_m));
		x_f = local_zeros(local_size(rhs_f));
		lagr  = local_zeros(local_size(rhs_f));

		UVector sol = blocks(x_m, x_f, lagr);

		Factorization<USparseMatrix, UVector> op;
		op.update(make_ref(A));
		bool ok = op.apply(rhs, sol);

		undo_blocks(sol, x_m, x_f, lagr);
		return ok;
	}

	bool SteadyFractureFlowSimulation::solve_staggered()
	{
		auto &V_m = matrix->space.subspace(0);
		auto &V_f = fracture_newtork->space.subspace(0);

	    Factorization<USparseMatrix, UVector> op_m;
	    op_m.update(make_ref(A_m));


	    assert(!empty(B));
	    assert(!empty(D));

	    L2TransferOperator t(make_ref(B), make_ref(D));

	    if(empty(x_f)) {
	        x_f = local_zeros(local_size(rhs_f));
	    }

	    UVector lagr_m = local_zeros(local_size(rhs_m));
	    UVector lagr_s = local_zeros(local_size(rhs_f));

	    UVector rhs_lagr_m, rhs_lagr_s;
	    UVector delta_lagr = local_zeros(local_size(rhs_f));

	    double dumping = 1.;

	    lagr_m.set(0.);

	    bool solved = false;
	    for(int i = 0; i < 20; ++i) {
	        apply_zero_boundary_conditions(V_m.dof_map(), lagr_m);
	        rhs_lagr_m = rhs_m + lagr_m;

	        op_m.apply(rhs_lagr_m, x_m);

	        x_f.set(0.);
	        t.apply(x_m, x_f);

	        apply_boundary_conditions(V_f, A_f, x_f);
	        lagr_s = rhs_f - A_f * x_f;
	        apply_zero_boundary_conditions(V_f.dof_map(), lagr_s);

	        double n_lagr_s = norm2(lagr_s);

	        disp(n_lagr_s);

	        if(n_lagr_s < 1e-14) {
	            solved = true;
	        }

	        delta_lagr.set(0);
	        t.apply_transpose(lagr_s, delta_lagr);
	        lagr_m += dumping * delta_lagr;
	    }

	    lagr = lagr_m;

	    return solved;
	}

	bool SteadyFractureFlowSimulation::solve_separate()
	{
		using SolverT = Factorization<USparseMatrix, UVector>;

	    SolverT().solve(A_m, rhs_m, x_m);
	    SolverT().solve(A_f, rhs_f, x_f);

	    assert(!empty(B));
	    assert(!empty(D));


	    UVector sol_transfered = local_zeros(local_size(D).get(0));
	    L2TransferOperator(make_ref(B), make_ref(D)).apply(x_m, sol_transfered);
	    lagr = sol_transfered - x_f;
	    return true;
	}

	void SteadyFractureFlowSimulation::write_output()
	{
		auto &V_m = matrix->space.subspace(0);
		auto &V_f = fracture_newtork->space.subspace(0);

		libMesh::Nemesis_IO io_m(matrix->mesh.mesh());
		libMesh::Nemesis_IO io_f(fracture_newtork->mesh.mesh());

		utopia::convert(x_m, *V_m.equation_system().solution);
		V_m.equation_system().solution->close();
		io_m.write_timestep(V_m.equation_system().name() + ".e", V_m.equation_systems(), 1, 0);

		utopia::convert(x_f, *V_f.equation_system().solution);
		V_f.equation_system().solution->close();
		io_f.write_timestep(V_f.equation_system().name() + ".e", V_f.equation_systems(), 1, 0);

		if(!lagrange_multiplier->empty()) {
		    libMesh::Nemesis_IO io_multiplier(lagrange_multiplier->mesh.mesh());
		    auto &L = lagrange_multiplier->space.subspace(0);
		    utopia::convert(lagr, *L.equation_system().solution);
		    L.equation_system().solution->close();
		    io_multiplier.write_timestep(L.equation_system().name() + ".e", L.equation_systems(), 1, 0);
		}
	}

	bool SteadyFractureFlowSimulation::run()
	{
		Chrono c;
		c.start();

		assemble_systems();
		init_coupling();

		bool ok = false;
		if(solve()) {
			write_output();
			ok = true;
		}

		c.stop();
		std::cout << "Overall time: " << c << std::endl;
		return ok;
	}

}
