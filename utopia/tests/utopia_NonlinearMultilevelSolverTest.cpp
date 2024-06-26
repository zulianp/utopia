#include "utopia.hpp"
#include "utopia_TestProblems.hpp"
#include "utopia_Testing.hpp"

namespace utopia {

#ifdef UTOPIA_ENABLE_PETSC
    class NonlinearBratuSolverTest {
    public:
        using SizeType = typename utopia::Traits<PetscVector>::SizeType;
        using Scalar = typename utopia::Traits<PetscVector>::Scalar;
        using Comm = typename utopia::Traits<PetscVector>::Communicator;

        explicit NonlinearBratuSolverTest(const SizeType& n_levels, bool remove_BC_contributions, bool verbose)
            : verbose_(verbose) {
            problem = MultiLevelTestProblem1D<PetscMatrix, PetscVector, Bratu1D<PetscMatrix, PetscVector>>(
                n_levels, 10, remove_BC_contributions);
        }

        void run() {
            UTOPIA_RUN_TEST(TR_test);
            UTOPIA_RUN_TEST(TR_constraint_test);

            UTOPIA_RUN_TEST(newton_MG_test);
            // UTOPIA_RUN_TEST(FAS_test);

            // UTOPIA_RUN_TEST(RMTR_test);
            // UTOPIA_RUN_TEST(RMTR_inf_test);
            // UTOPIA_RUN_TEST(RMTR_inf_bound_test);
        }

        void TR_test() {
            Bratu1D<PetscMatrix, PetscVector> fun(1000);
            PetscVector x = fun.initial_guess();

            auto subproblem = std::make_shared<utopia::SteihaugToint<PetscMatrix, PetscVector, HOMEMADE>>();
            subproblem->set_preconditioner(std::make_shared<InvDiagPreconditioner<PetscMatrix, PetscVector>>());
            subproblem->atol(1e-14);
            subproblem->max_it(1000);

            TrustRegion<PetscMatrix, PetscVector> tr_solver(subproblem);
            tr_solver.atol(1e-9);
            tr_solver.rtol(1e-9);
            tr_solver.stol(1e-10);
            tr_solver.max_it(500);
            tr_solver.verbose(verbose_);
            tr_solver.solve(fun, x);
            // disp(x);
        }

        void TR_constraint_test() {
            Bratu1D<PetscMatrix, PetscVector> fun(problem.n_dofs(1));
            PetscVector x = fun.initial_guess();
            // fun.apply_bc_to_initial_guess(x);

            PetscVector ub, lb;
            fun.generate_constraints(lb, ub);
            auto box = make_box_constaints(make_ref(lb), make_ref(ub));

            auto lsolver = std::make_shared<LUDecomposition<PetscMatrix, PetscVector>>();
            auto qp_solver = std::make_shared<utopia::TaoQPSolver<PetscMatrix, PetscVector>>(lsolver);
            qp_solver->atol(1e-11);

            TrustRegionVariableBound<PetscMatrix, PetscVector> tr_solver(qp_solver);
            tr_solver.set_box_constraints(box);
            tr_solver.atol(1e-8);
            tr_solver.rtol(1e-10);
            tr_solver.stol(1e-10);
            tr_solver.verbose(verbose_);

            tr_solver.solve(fun, x);

            PetscVector x2 = fun.initial_guess();
            fun.apply_bc_to_initial_guess(x2);
            auto qp_solver2 = std::make_shared<utopia::SemismoothNewton<PetscMatrix, PetscVector>>(lsolver);
            qp_solver2->atol(1e-11);
            tr_solver.set_trust_region_strategy(qp_solver2);
            tr_solver.solve(fun, x2);

            // takes forever
            // PetscVector x3 = fun.initial_guess();
            // fun.apply_bc_to_initial_guess(x3);
            // auto qp_solver3 =  std::make_shared<utopia::ProjectedGradient<PetscMatrix, PetscVector> >();
            // qp_solver3->atol(1e-11);
            // tr_solver.set_trust_region_strategy(qp_solver3);
            // tr_solver.solve(fun, x3);

            PetscVector x4 = fun.initial_guess();
            fun.apply_bc_to_initial_guess(x4);
            auto qp_solver4 = std::make_shared<utopia::ProjectedConjugateGradient<PetscMatrix, PetscVector>>();
            qp_solver4->atol(1e-11);
            tr_solver.set_trust_region_strategy(qp_solver4);
            tr_solver.solve(fun, x4);

            utopia_test_assert(approxeq(x, x2));
            utopia_test_assert(approxeq(x, x4));
        }

        void newton_MG_test() {
            auto lsolver = std::make_shared<utopia::BiCGStab<PetscMatrix, PetscVector>>();
            Newton<utopia::PetscMatrix, utopia::PetscVector> newton(lsolver);

            auto direct_solver = std::make_shared<LUDecomposition<PetscMatrix, PetscVector>>();
            auto gs = std::make_shared<GaussSeidel<PetscMatrix, PetscVector>>();
            auto multigrid = std::make_shared<Multigrid<PetscMatrix, PetscVector>>(gs, direct_solver);

            auto funs = problem.get_functions();

            PetscVector x_0;
            funs[problem.n_levels() - 1]->get_eq_constrains_values(x_0);

            multigrid->set_transfer_operators(problem.get_transfer());
            multigrid->must_generate_masks(false);
            multigrid->fix_semidefinite_operators(true);
            multigrid->verbose(false);
            multigrid->atol(1e-11);
            multigrid->rtol(1e-11);

            newton.set_linear_solver(multigrid);
            newton.verbose(verbose_);
            newton.atol(1e-9);
            newton.rtol(1e-10);
            newton.solve(*funs.back(), x_0);
        }

        // void FAS_test() {
        //     // intial guess
        //     Comm comm;
        //     PetscVector x;
        //     x.values(comm, problem.n_dofs(problem.n_levels() - 1), 0.0);

        //     std::vector<std::shared_ptr<ExtendedFunction<PetscMatrix, PetscVector>>>
        //     level_functions(problem.n_levels);

        //     for (auto l = 0; l < problem.n_levels; l++) {
        //         auto fun = std::make_shared<Bratu1D<PetscMatrix, PetscVector>>(problem.n_dofs(l));
        //         level_functions[l] = fun;

        //         // making sure that fine level IG is feasible
        //         if (l + 1 == problem.n_levels) fun->apply_bc_to_initial_guess(x);
        //     }

        //     auto direct_solver = std::make_shared<LUDecomposition<PetscMatrix, PetscVector>>();
        //     auto coarse_solver = std::make_shared<Newton<utopia::PetscMatrix, utopia::PetscVector>>(direct_solver);
        //     auto strategy = std::make_shared<utopia::Backtracking<utopia::PetscVector>>();
        //     coarse_solver->set_line_search_strategy(strategy);
        //     coarse_solver->atol(1e-7);
        //     coarse_solver->max_it(1);

        //     // subject to change
        //     auto smoother = std::make_shared<NonLinearJacobi<PetscMatrix, PetscVector>>();
        //     smoother->relaxation_parameter(0.3);
        //     // auto smoother = std::make_shared<NonLinearGMRES<PetscMatrix, PetscVector> >();

        //     auto fas = std::make_shared<FAS<PetscMatrix, PetscVector>>(problem.n_levels, smoother, coarse_solver);
        //     fas->set_transfer_operators(problem.prolongations, problem.restrictions, problem.restrictions);

        //     fas->pre_smoothing_steps(3);
        //     fas->post_smoothing_steps(3);
        //     fas->verbose(problem.verbose);
        //     fas->atol(1e-8);
        //     fas->rtol(1e-10);
        //     fas->max_it(100);

        //     fas->set_functions(level_functions);
        //     fas->solve(x);
        // }

        // void RMTR_test() {
        //     // intial guess
        //     PetscVector x = values(problem.n_dofs(problem.n_levels() - 1), 0.0);

        //     std::vector<std::shared_ptr<ExtendedFunction<PetscMatrix, PetscVector>>>
        //     level_functions(problem.n_levels);

        //     for (auto l = 0; l < problem.n_levels; l++) {
        //         auto fun = std::make_shared<Bratu1D<PetscMatrix, PetscVector>>(problem.n_dofs(l));
        //         level_functions[l] = fun;

        //         // making sure that fine level IG is feasible
        //         if (l + 1 == problem.n_levels()) fun->apply_bc_to_initial_guess(x);
        //     }

        //     auto tr_strategy_coarse = std::make_shared<utopia::SteihaugToint<PetscMatrix, PetscVector, HOMEMADE>>();
        //     tr_strategy_coarse->atol(1e-9);
        //     // tr_strategy_coarse->verbose(true);
        //     auto tr_strategy_fine = std::make_shared<utopia::SteihaugToint<PetscMatrix, PetscVector, HOMEMADE>>();
        //     tr_strategy_fine->atol(1e-9);
        //     // tr_strategy_fine->verbose(true);

        //     // we should apply BC conditions in symmetric way... ASAP...
        //     tr_strategy_coarse->set_preconditioner(std::make_shared<IdentityPreconditioner<PetscVector>>());
        //     // tr_strategy_fine->set_preconditioner(std::make_shared<InvDiagPreconditioner<PetscMatrix, PetscVector>
        //     // >());

        //     tr_strategy_coarse->set_preconditioner(std::make_shared<IdentityPreconditioner<PetscVector>>());
        //     // tr_strategy_fine->set_preconditioner(std::make_shared<InvDiagPreconditioner<PetscMatrix, PetscVector>
        //     // >());

        //     // auto rmtr = std::make_shared<RMTR<PetscMatrix, PetscVector, SECOND_ORDER>  >(tr_strategy_coarse,
        //     // tr_strategy_fine);
        //     auto rmtr = std::make_shared<RMTR<PetscMatrix, PetscVector, FIRST_ORDER>>(problem.n_levels);
        //     rmtr->set_transfer_operators(problem.prolongations(), problem.restrictions());

        //     rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
        //     rmtr->set_fine_tr_strategy(tr_strategy_fine);

        //     rmtr->max_it(50);
        //     rmtr->max_coarse_it(1);
        //     rmtr->max_sucessful_smoothing_it(3);

        //     rmtr->delta0(10);
        //     rmtr->atol(1e-5);
        //     rmtr->rtol(1e-10);
        //     rmtr->set_grad_smoothess_termination(0.000001);

        //     rmtr->verbose(problem.verbose);
        //     rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);
        //     rmtr->set_functions(level_functions);

        //     rmtr->solve(x);
        // }

        // void RMTR_inf_test() {
        //     // intial guess
        //     PetscVector x = values(problem.n_dofs(problem.n_levels() - 1), 0.0);

        //     std::vector<std::shared_ptr<ExtendedFunction<PetscMatrix, PetscVector>>>
        //     level_functions(problem.n_levels); for (auto l = 0; l < problem.n_levels; l++) {
        //         auto fun = std::make_shared<Bratu1D<PetscMatrix, PetscVector>>(problem.n_dofs(l));
        //         level_functions[l] = fun;

        //         // making sure that fine level IG is feasible
        //         if (l + 1 == problem.n_levels) fun->apply_bc_to_initial_guess(x);
        //     }

        //     auto lsolver = std::make_shared<LUDecomposition<PetscMatrix, PetscVector>>();
        //     auto tr_strategy_fine = std::make_shared<utopia::TaoQPSolver<PetscMatrix, PetscVector>>(lsolver);
        //     tr_strategy_fine->set_linear_solver(std::make_shared<GMRES<PetscMatrix, PetscVector>>(PCJACOBI));

        //     auto tr_strategy_coarse = std::make_shared<utopia::TaoQPSolver<PetscMatrix, PetscVector>>(lsolver);
        //     tr_strategy_coarse->set_linear_solver(std::make_shared<GMRES<PetscMatrix, PetscVector>>(PCLU));

        //     auto rmtr = std::make_shared<RMTR_inf<PetscMatrix, PetscVector, SECOND_ORDER>>(problem.n_levels);
        //     rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
        //     rmtr->set_fine_tr_strategy(tr_strategy_fine);

        //     rmtr->set_transfer_operators(problem.prolongations, problem.restrictions);

        //     rmtr->max_it(1000);
        //     rmtr->max_coarse_it(1);
        //     rmtr->max_sucessful_smoothing_it(1);
        //     rmtr->delta0(1);
        //     rmtr->atol(1e-5);
        //     rmtr->rtol(1e-10);
        //     rmtr->set_grad_smoothess_termination(0.000001);

        //     rmtr->verbose(problem.verbose);
        //     // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
        //     rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);

        //     rmtr->set_functions(level_functions);
        //     rmtr->solve(x);
        // }

        // void RMTR_inf_bound_test() {
        //     // intial guess
        //     PetscVector x = values(problem.n_dofs[problem.n_levels() - 1], 0.0);

        //     // upper, lower bound...
        //     PetscVector ub, lb;
        //     std::vector<std::shared_ptr<ExtendedFunction<PetscMatrix, PetscVector>>>
        //     level_functions(problem.n_levels); for (auto l = 0; l < problem.n_levels; l++) {
        //         auto fun = std::make_shared<Bratu1D<PetscMatrix, PetscVector>>(problem.n_dofs(l));
        //         level_functions[l] = fun;

        //         // making sure that fine level IG is feasible
        //         if (l + 1 == problem.n_levels) {
        //             fun->apply_bc_to_initial_guess(x);
        //             fun->generate_constraints(lb, ub, -10, 0.1);
        //         }
        //     }

        //     // Utopia::instance().set("log_output_path", "benchmark.csv");

        //     auto lsolver = std::make_shared<LUDecomposition<PetscMatrix, PetscVector>>();
        //     auto tr_strategy_fine = std::make_shared<utopia::TaoQPSolver<PetscMatrix, PetscVector>>(lsolver);
        //     tr_strategy_fine->set_linear_solver(std::make_shared<GMRES<PetscMatrix, PetscVector>>(PCJACOBI));
        //     tr_strategy_fine->verbose(false);

        //     auto tr_strategy_coarse = std::make_shared<utopia::TaoQPSolver<PetscMatrix, PetscVector>>(lsolver);
        //     tr_strategy_coarse->set_linear_solver(std::make_shared<GMRES<PetscMatrix, PetscVector>>(PCLU));
        //     // tr_strategy_coarse->verbose(true);
        //     tr_strategy_coarse->verbose(false);

        //     auto rmtr = std::make_shared<RMTR_inf<PetscMatrix, PetscVector, SECOND_ORDER>>(problem.n_levels);
        //     rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
        //     rmtr->set_fine_tr_strategy(tr_strategy_fine);

        //     rmtr->set_transfer_operators(problem.prolongations, problem.restrictions);

        //     rmtr->max_it(1000);
        //     rmtr->max_coarse_it(1);
        //     rmtr->max_sucessful_smoothing_it(1);
        //     rmtr->delta0(1);
        //     rmtr->atol(1e-5);
        //     rmtr->rtol(1e-10);
        //     rmtr->set_grad_smoothess_termination(0.000001);

        //     rmtr->verbose(problem.verbose);
        //     // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
        //     rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);

        //     rmtr->set_functions(level_functions);

        //     auto box = make_box_constaints(make_ref(lb), make_ref(ub));
        //     rmtr->set_box_constraints(box);
        //     rmtr->solve(x);
        // }

    private:
        MultilevelTestProblemBase<PetscMatrix, PetscVector> problem;
        bool verbose_;
    };

#endif  // UTOPIA_ENABLE_PETSC

    static void non_linear_multilevel() {
#ifdef UTOPIA_ENABLE_PETSC
        NonlinearBratuSolverTest(3, true, false).run();
#endif
    }

    UTOPIA_REGISTER_TEST_FUNCTION(non_linear_multilevel);
}  // namespace utopia
