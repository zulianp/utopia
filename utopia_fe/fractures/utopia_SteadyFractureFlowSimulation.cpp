#include "utopia_SteadyFractureFlowSimulation.hpp"

#include "utopia.hpp"
#include "utopia_Socket.hpp"
#include "utopia_FractureFlowUtils.hpp"
#include "utopia_TransferUtils.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_SPStaticCondensation.hpp"
#include "utopia_SPStaticCondensationKrylov.hpp"


#include <iostream>


namespace utopia {

    void SteadyFractureFlowSimulation::describe(std::ostream &os) const
    {
        matrix->describe(os);
        fracture_network->describe(os);
        if(lagrange_multiplier) lagrange_multiplier->describe(os);

        std::cout << "solve_strategy                " << solve_strategy << std::endl;
        std::cout << "use_mg                        " << use_mg << std::endl;
        std::cout << "mg_sweeps                     " << mg_sweeps << std::endl;
        std::cout << "mg_levels                     " << mg_levels << std::endl;
        std::cout << "plot_matrix                   " << plot_matrix << std::endl;
        std::cout << "write_operators_to_disk       " << write_operators_to_disk << std::endl;
        std::cout << "normal_hydraulic_conductivity	" << normal_hydraulic_conductivity << std::endl;
    }

    SteadyFractureFlowSimulation::SteadyFractureFlowSimulation(libMesh::Parallel::Communicator &comm)
    {
        matrix 				= utopia::make_unique<FractureFlow>(comm);
        fracture_network 	= utopia::make_unique<FractureFlow>(comm);
        lagrange_multiplier = utopia::make_unique<FractureFlow>(comm);

        solve_strategy = "monolithic";
        use_mg = false;
        mg_sweeps = 1;
        mg_levels = 5;
        plot_matrix = false;
        write_operators_to_disk = false;
        normal_hydraulic_conductivity = 1.;
        use_interpolation = false;
        use_biorth = false;
    }

    void SteadyFractureFlowSimulation::read(Input &is)
    {
        Chrono c;
        c.start();

        is.get("master", 	 *matrix);
        is.get("slave",  	 *fracture_network);
        is.get("multiplier", *lagrange_multiplier);

        is.get("solve-strategy", solve_strategy);
        is.get("use-mg", use_mg);
        is.get("mg-sweeps", mg_sweeps);
        is.get("mg-levels", mg_levels);
        is.get("plot-matrix", plot_matrix);
        is.get("normal-hydraulic-conductivity", normal_hydraulic_conductivity);
        is.get("write-operators-to-disk", write_operators_to_disk);
        is.get("use-interpolation", use_interpolation);
        is.get("use-biorth", use_biorth);

        auto subdomain_fun = utopia::make_unique<UISubdomainFunction<double>>();

        is.get("normal-hydraulic-conductivity-blocks", *subdomain_fun);

        if(!subdomain_fun->good()) {
            normal_hydraulic_conductivity_blocks = std::make_shared<UIConstantFunction<double>>(1.);
        } else {

            if(!subdomain_fun->has_default()) {
                subdomain_fun->set_default(utopia::make_unique<UIConstantFunction<double>>(1.));
            }

            normal_hydraulic_conductivity_blocks = std::move(subdomain_fun);
        }

        if(plot_matrix) {
            plot_mesh(matrix->mesh.mesh(), "matrix");
        }

        this->describe(std::cout);

        auto &V_m = matrix->space.subspace(0);
        auto &V_f = fracture_network->space.subspace(0);

        auto aux_matrix = utopia::make_unique<FractureFlowAuxSystem>(V_m);
        aux_matrix->sample(matrix->sampler);

        auto aux_fracture_network = utopia::make_unique<FractureFlowAuxSystem>(V_f);
        aux_fracture_network->sample(fracture_network->sampler);

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
        auto &V_f = fracture_network->space.subspace(0);

        auto u_m = trial(V_m);
        auto v_m = test(V_m);

        auto u_s = trial(V_f);
        auto v_f = test(V_f);

        auto eq_m = inner(matrix->diffusion_tensor * grad(u_m), ctx_fun(matrix->sampler) * grad(v_m)) * dX;
        auto eq_f = inner(fracture_network->diffusion_tensor  * grad(u_s), ctx_fun(fracture_network->sampler)  * grad(v_f)) * dX;

        // auto eq_m = inner(grad(u_m), ctx_fun(matrix->sampler) * grad(v_m)) * dX;
        // auto eq_f = inner(grad(u_s), ctx_fun(fracture_network->sampler)  * grad(v_f)) * dX;


        utopia::assemble(eq_m, A_m);
        utopia::assemble(eq_f, A_f);

        x_m = local_zeros(V_m.dof_map().n_local_dofs());
        x_f = local_zeros(V_f.dof_map().n_local_dofs());

        matrix->forcing_function->eval(x_m, rhs_m);
        fracture_network->forcing_function->eval(x_f, rhs_f);

        double norm_f_f = norm2(rhs_f);
        double norm_f_m = norm2(rhs_m);

        std::cout << "forcing_functions: " << norm_f_m << " " << norm_f_f << std::endl;

        matrix->apply_weak_BC(A_m, rhs_m);
        fracture_network->apply_weak_BC(A_f, rhs_f);

        if(write_operators_to_disk) {
            write("A_m_neu.m", A_m);
            write("A_f_neu.m", A_f);
        }

        apply_boundary_conditions(V_m.dof_map(), A_m, rhs_m);
        apply_boundary_conditions(V_f.dof_map(), A_f, rhs_f);


        c.stop();
        std::cout << "model assemly time: " << c << std::endl;

        if(write_operators_to_disk) {
            write("A_m.m", A_m);
            write("A_f.m", A_f);
        }
    }

    template<class Space>
    static void project_function(const std::shared_ptr<UIFunction<double>> &f, Space &V, UVector &vec)
    {
        UVector lumped_mass;
        UVector vec_h;
        auto v = test(V);
        utopia::assemble(inner(coeff(1.), v) * dX, lumped_mass);
        utopia::assemble(inner(ctx_fun(f), v) * dX, vec_h);
        vec = e_mul(vec_h, 1./lumped_mass);
    }

    void SteadyFractureFlowSimulation::init_coupling()
    {
        Chrono c;
        c.start();

        auto &V_m = matrix->space.subspace(0);
        auto &V_f = fracture_network->space.subspace(0);

        if(lagrange_multiplier->empty()) {
            if(use_interpolation) {
                assemble_interpolation(V_m, V_f, B, D);
            } else {
                assemble_projection(V_m, V_f, B, D, use_biorth);
            }
        } else {
            auto &V_l = lagrange_multiplier->space.subspace(0);
            assemble_projection(V_m, V_f, V_l, B, D);

        }

        D *= -1.;

        D_t = transpose(D);
        B_t = transpose(B);

        set_zero_at_constraint_rows(V_m.dof_map(), B_t);
        set_zero_at_constraint_rows(V_f.dof_map(), D_t);



        auto constant_function = std::dynamic_pointer_cast<UIConstantFunction<double>>(normal_hydraulic_conductivity_blocks);
        if(constant_function) {
            kappa_B_t = B_t;
            kappa_B_t *= (normal_hydraulic_conductivity * constant_function->value());

            std::cout << "normal-hydraulic-conductivity-blocks = " << constant_function->value() << std::endl;
        } else {
            //perform l2-projection of kappa

            //either on lagrange mult space or slave space
            if(lagrange_multiplier->empty()) {
                UVector kappa;
                project_function(normal_hydraulic_conductivity_blocks, V_f, kappa);
                USparseMatrix K = diag(kappa);
                kappa_B_t = B_t * K;
                kappa_B_t *= normal_hydraulic_conductivity;

                double min_kappa = min(kappa), max_kappa = max(kappa);
                std::cout << "normal-hydraulic-conductivity-blocks  in [" << min_kappa << ", " << max_kappa << "]" << std::endl;
            } else {
                auto &V_l = lagrange_multiplier->space.subspace(0);
                UVector kappa;
                project_function(normal_hydraulic_conductivity_blocks, V_l, kappa);
                USparseMatrix K = diag(kappa);
                kappa_B_t = B_t * K;
                kappa_B_t *= normal_hydraulic_conductivity;

                double min_kappa = min(kappa), max_kappa = max(kappa);
                std::cout << "normal-hydraulic-conductivity-blocks  in [" << min_kappa << ", " << max_kappa << "]" << std::endl;
            }

        }

        if(write_operators_to_disk) {
            write("B.m", B);
            write("D.m", D);
        }

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
        } else if(solve_strategy == "static-condensation") {
            ok = solve_monolithic_static_condenstation();
        } else {
            ok = solve_monolithic();
        }

        matrix->post_process(x_m);
        fracture_network->post_process(x_f);

        c.stop();
        std::cout << "Solver time: " << c << std::endl;
        return ok;
    }

    /*
        A_m 0   B_t
        0	A_f kappa * D_t
        B   kappa * D   0

    */

    bool SteadyFractureFlowSimulation::solve_cg_dual()
    {
        auto &V_m = matrix->space.subspace(0);
        auto &V_f = fracture_network->space.subspace(0);

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
            make_ref(kappa_B_t),
            make_ref(D_t)
        );

        return solver.apply(rhs_m, rhs_f, x_m, x_f, lagr);
    }

    bool SteadyFractureFlowSimulation::solve_monolithic()
    {
        USparseMatrix A = Blocks<USparseMatrix>(3, 3,
        {
            make_ref(A_m), nullptr, make_ref(kappa_B_t),
            nullptr, make_ref(A_f), make_ref(D_t),
            make_ref(B), make_ref(D), nullptr
        });

        UVector z = local_zeros(local_size(B).get(0));
        UVector rhs = blocks(rhs_m, rhs_f, z);

        x_m  = local_zeros(local_size(rhs_m));
        x_f  = local_zeros(local_size(rhs_f));
        lagr = local_zeros(local_size(z));

        UVector sol = blocks(x_m, x_f, lagr);


        if(write_operators_to_disk) {
            write("S_moonolithic.m", A);
        }

        Factorization<USparseMatrix, UVector> op;
        op.update(make_ref(A));
        op.describe(std::cout);

        bool ok = op.apply(rhs, sol);

        undo_blocks(sol, x_m, x_f, lagr);
        return ok;
    }

    // bool SteadyFractureFlowSimulation::solve_monolithic_static_condenstation()
    // {
    // 	// std::cout << "solve_monolithic_static_condenstation" << std::endl;
    // 	auto &V_m = matrix->space.subspace(0);
    // 	auto &V_f = fracture_network->space.subspace(0);

    // 	//Remove the - from D
    // 	UVector d_inv = (-1.)/sum(D, 1);
    // 	USparseMatrix D_inv = diag(d_inv);
    // 	T = D_inv * B;

    // 	apply_boundary_conditions(V_m.dof_map(), A_m, rhs_m);
    // 	apply_boundary_conditions(V_f.dof_map(), A_f, rhs_f);

    // 	USparseMatrix S = A_m + transpose(T) * A_f * T;
    // 	UVector rhs = rhs_m + transpose(T) * rhs_f;
    // 	// apply_boundary_conditions(V_m.dof_map(), S, rhs);

    //     if(use_mg) {
    //     	auto mg = make_mg_solver(V_m, mg_levels);
    //     	if(!mg->solve(S, rhs, x_m)) {
    //     		return false;
    //     	}

    //     } else {
    // 		Factorization<USparseMatrix, UVector> op_m;
    // 		if(!op_m.solve(S, rhs, x_m)) {
    // 			return false;
    // 		}
    // 	}

    // 	x_f = T * x_m;
    // 	lagr = D_inv * (A_f * (T * x_m) - rhs_f);
    // 	return true;
    // }

    bool SteadyFractureFlowSimulation::solve_monolithic_static_condenstation()
    {
        auto &V_m = matrix->space.subspace(0);
        auto &V_f = fracture_network->space.subspace(0);

        apply_boundary_conditions(V_m.dof_map(), A_m, rhs_m);
        apply_boundary_conditions(V_f.dof_map(), A_f, rhs_f);

        //also add dual-basis case
        if(use_interpolation || use_biorth) {
            //Remove the - from D
            UVector d_inv = (-1.)/sum(D, 1);
            USparseMatrix D_inv = diag(d_inv);
            T = D_inv * B;

            SPStaticCondensation<USparseMatrix, UVector> sp(std::make_shared<Factorization<USparseMatrix, UVector>>());

            if(use_mg) {
                sp.linear_solver(make_mg_solver(V_m, mg_levels));
            }

            sp.update(make_ref(A_m), make_ref(A_f), make_ref(T));

            // apply_boundary_conditions(V_m.dof_map(), sp.get_operator(), x_m);
            
            if(!sp.apply(rhs_m, rhs_f, x_m, x_f)) {
                return false;
            }

            if(write_operators_to_disk) {
                write("S_static.m", sp.get_operator());
            }

            lagr = D_inv * (A_f * (T * x_m) - rhs_f);
            return true;
        } else {

            SPStaticCondensationKrylov<USparseMatrix, UVector> sp;

            // auto cg = std::make_shared<BiCGStab<USparseMatrix, UVector, HOMEMADE>>();
            auto mf_solver = std::make_shared<ProjectedGradient<USparseMatrix, UVector, HOMEMADE>>();
            mf_solver->verbose(true);
            mf_solver->max_it(60000);

            sp.linear_solver(mf_solver);
            sp.coupling_op_solver(std::make_shared<Factorization<USparseMatrix, UVector>>());


            // if(use_mg) {
            // 	sp.linear_solver(make_mg_solver(V_m, mg_levels));
            // }
            USparseMatrix m_D = D;
            m_D *= -1.;

            sp.update(make_ref(A_m), make_ref(A_f), make_ref(B), make_ref(m_D));

            if(!sp.apply(rhs_m, rhs_f, x_m, x_f)) {
                return false;
            }

            UVector temp;
            Factorization<USparseMatrix, UVector> fact;
            bool ok = fact.solve(m_D, B * x_m, temp);
            temp = A_f * temp - rhs_f;
            return fact.apply(temp, lagr) && ok;
        }
    }

    bool SteadyFractureFlowSimulation::solve_staggered()
    {
        auto &V_m = matrix->space.subspace(0);
        auto &V_f = fracture_network->space.subspace(0);

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
        L2TransferOperator t_op(make_ref(B), make_ref(D));
        t_op.init();
        t_op.apply(x_m, sol_transfered);
        lagr = sol_transfered - x_f;
        return true;
    }

    template<class IO>
    void SteadyFractureFlowSimulation::write_output_generic()
    {
        auto &V_m = matrix->space.subspace(0);
        auto &V_f = fracture_network->space.subspace(0);

        IO io_m(matrix->mesh.mesh());
        IO io_f(fracture_network->mesh.mesh());

        utopia::convert(x_m, *V_m.equation_system().solution);
        V_m.equation_system().solution->close();
        io_m.write_equation_systems(V_m.equation_system().name() + ".e", V_m.equation_systems());

        utopia::convert(x_f, *V_f.equation_system().solution);
        V_f.equation_system().solution->close();
        io_f.write_equation_systems(V_f.equation_system().name() + ".e", V_f.equation_systems());

        if(!lagrange_multiplier->empty()) {
            IO io_multiplier(lagrange_multiplier->mesh.mesh());
            auto &L = lagrange_multiplier->space.subspace(0);
            utopia::convert(lagr, *L.equation_system().solution);
            L.equation_system().solution->close();
            io_multiplier.write_equation_systems(L.equation_system().name() + ".e", L.equation_systems());
        }
    }

    void SteadyFractureFlowSimulation::write_output()
    {
        auto &V_m = matrix->space.subspace(0);
        auto &V_f = fracture_network->space.subspace(0);

        if(V_m.mesh().comm().size() == 1) {
            write_output_generic<libMesh::ExodusII_IO>();
        } else {
            write_output_generic<libMesh::Nemesis_IO>();
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
