#include "utopia.hpp"
#include "utopia_SolverTest.hpp"
#include "test_problems/utopia_TestProblems.hpp"
#include "test_problems/utopia_BratuMultilevelTestProblem.hpp"

namespace utopia
{


#ifdef  WITH_PETSC
    class NonlinearBratuSolverTest {
    public:

        typedef UTOPIA_SIZE_TYPE(DVectord) SizeType;
        typedef UTOPIA_SCALAR(DVectord) Scalar;

        NonlinearBratuSolverTest(const SizeType & n_levels = 2, bool remove_BC_contributions = false, bool verbose = false):
        problem(n_levels, remove_BC_contributions, verbose)
        {

        }

        void run()
        {
            // UTOPIA_RUN_TEST(TR_test);
            // UTOPIA_RUN_TEST(TR_constraint_test);

            // UTOPIA_RUN_TEST(newton_MG_test);
            // UTOPIA_RUN_TEST(FAS_test);

            // UTOPIA_RUN_TEST(RMTR_test);
            // UTOPIA_RUN_TEST(RMTR_inf_test);
            // UTOPIA_RUN_TEST(RMTR_inf_bound_test);


            UTOPIA_RUN_TEST(PETSC_test);

        }

        void PETSC_test()
        {
            std::cout<<"----- petsc test.... \n"; 

            auto n = 150; 

            // AppCtxBratu2D app_data; 
            // app_data.lambda  = 1.0; 
            // app_data.mms_solution = Bratu2DMMSSolution; 
            // app_data.mms_forcing = Bratu2DMMSForcing;            


            // DM da_coarse, da_fine; 
            // DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,DMDA_STENCIL_STAR, n, n,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL, &da_coarse);
            // DMSetUp(da_coarse);
            // DMDASetUniformCoordinates(da_coarse, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
            // DMRefine(da_coarse,PETSC_COMM_WORLD, &da_fine); 
            // DMDestroy(&da_coarse); 


            // Bratu2D<DSMatrixd, DVectord> fun(da_fine, app_data);
            Bratu2D<DSMatrixd, DVectord> fun(n, 1.0);
            DVectord x = fun.initial_guess(); 
            fun.describe(); 
            

            auto subproblem = std::make_shared<utopia::Lanczos<DSMatrixd, DVectord> >();
            subproblem->pc_type("bjacobi"); 
            subproblem->atol(1e-14);
            subproblem->max_it(100000);
            TrustRegion<DSMatrixd, DVectord> tr_solver(subproblem);
            tr_solver.atol(1e-5);
            tr_solver.rtol(1e-10);
            tr_solver.stol(1e-10);
            tr_solver.verbose(true);
            
            // x = 0.0*x; 
            tr_solver.solve(fun, x);

            fun.output_to_VTK(x);
        }




        void TR_test()
        {
            Bratu1D<DSMatrixd, DVectord> fun(problem.n_coarse);
            DVectord x = values(problem.n_coarse, 1.0);
            fun.apply_bc_to_initial_guess(x);

            auto subproblem = std::make_shared<utopia::SteihaugToint<DSMatrixd, DVectord> >();
            TrustRegion<DSMatrixd, DVectord> tr_solver(subproblem);
            tr_solver.atol(1e-10);
            tr_solver.rtol(1e-10);
            tr_solver.stol(1e-10);
            tr_solver.verbose(problem.verbose);
            tr_solver.solve(fun, x);
        }

        void TR_constraint_test()
        {
            Bratu1D<DSMatrixd, DVectord> fun(problem.n_coarse);
            DVectord x = values(problem.n_coarse, 1.0);
            fun.apply_bc_to_initial_guess(x);

            DVectord ub, lb;
            fun.generate_constraints(lb, ub);
            auto box = make_box_constaints(make_ref(lb), make_ref(ub));

            auto lsolver = std::make_shared<LUDecomposition<DSMatrixd, DVectord> >();
            auto qp_solver =  std::make_shared<utopia::TaoQPSolver<DSMatrixd, DVectord> >(lsolver);
            qp_solver->atol(1e-11);


            TrustRegionVariableBound<DSMatrixd, DVectord>  tr_solver(qp_solver);
            tr_solver.set_box_constraints(box);
            tr_solver.atol(1e-6);
            tr_solver.rtol(1e-10);
            tr_solver.stol(1e-10);
            tr_solver.verbose(problem.verbose);


            tr_solver.solve(fun, x);

            DVectord x2 = values(problem.n_coarse, 1.0);
            fun.apply_bc_to_initial_guess(x2);
            auto qp_solver2 =  std::make_shared<utopia::SemismoothNewton<DSMatrixd, DVectord> >(lsolver);
            qp_solver2->atol(1e-11);
            tr_solver.set_trust_region_strategy(qp_solver2);
            tr_solver.solve(fun, x2);

            DVectord x3 = values(problem.n_coarse, 1.0);
            fun.apply_bc_to_initial_guess(x3);
            auto qp_solver3 =  std::make_shared<utopia::ProjectedGradient<DSMatrixd, DVectord> >();
            qp_solver3->atol(1e-11);
            tr_solver.set_trust_region_strategy(qp_solver3);
            tr_solver.solve(fun, x3);

            DVectord x4 = values(problem.n_coarse, 1.0);
            fun.apply_bc_to_initial_guess(x4);
            auto qp_solver4 =  std::make_shared<utopia::ProjectedConjugateGradient<DSMatrixd, DVectord> >();
            qp_solver4->atol(1e-11);
            tr_solver.set_trust_region_strategy(qp_solver4);
            tr_solver.solve(fun, x4);

            utopia_test_assert(approxeq(x, x2));
            utopia_test_assert(approxeq(x3, x4));

        }

        void newton_MG_test()
        {
            Bratu1D<DSMatrixd, DVectord> fun(problem.n_dofs[problem.n_levels - 1]);
            DVectord x = values(problem.n_dofs[problem.n_levels - 1], 1.0);
            fun.apply_bc_to_initial_guess(x);

            auto lsolver = std::make_shared<utopia::BiCGStab<DSMatrixd, DVectord> >();
            Newton<utopia::DSMatrixd, utopia::DVectord> newton(lsolver);

            auto direct_solver = std::make_shared<LUDecomposition<DSMatrixd, DVectord>>();
            auto gs = std::make_shared<GaussSeidel<DSMatrixd, DVectord> >();
            auto multigrid = std::make_shared<Multigrid<DSMatrixd, DVectord>  >(gs, direct_solver);

            multigrid->set_transfer_operators(problem.prolongations);
            multigrid->must_generate_masks(false);
            multigrid->fix_semidefinite_operators(true);
            multigrid->verbose(false);
            multigrid->atol(1e-11);

            newton.set_linear_solver(multigrid);
            newton.verbose(problem.verbose);
            newton.atol(1e-9);
            newton.rtol(1e-10);
            newton.solve(fun, x);
        }


        void FAS_test()
        {
            // intial guess
            DVectord x = values(problem.n_dofs[problem.n_levels -1 ], 0.0);

            std::vector<std::shared_ptr<ExtendedFunction<DSMatrixd, DVectord> > >  level_functions(problem.n_levels);


            for(auto l=0; l < problem.n_levels; l++)
            {
                auto fun = std::make_shared<Bratu1D<DSMatrixd, DVectord> >(problem.n_dofs[l]);
                level_functions[l] = fun;

                // making sure that fine level IG is feasible
                if(l+1 == problem.n_levels)
                    fun->apply_bc_to_initial_guess(x);
            }


            auto direct_solver = std::make_shared<LUDecomposition<DSMatrixd, DVectord>>();
            auto coarse_solver = std::make_shared<Newton<utopia::DSMatrixd, utopia::DVectord> >(direct_solver);
               auto strategy = std::make_shared<utopia::Backtracking<utopia::DVectord>>();
            coarse_solver->set_line_search_strategy(strategy);
            coarse_solver->atol(1e-7);
            coarse_solver->max_it(1);

            // subject to change
            auto smoother = std::make_shared<NonLinearJacobi<DSMatrixd, DVectord> >();
            smoother->relaxation_parameter(0.3);
            // auto smoother = std::make_shared<NonLinearGMRES<DSMatrixd, DVectord> >();


            auto fas = std::make_shared<FAS<DSMatrixd, DVectord>  >(problem.n_levels, smoother, coarse_solver);
            fas->set_transfer_operators(problem.prolongations, problem.restrictions, problem.restrictions);

            fas->pre_smoothing_steps(3);
            fas->post_smoothing_steps(3);
            fas->verbose(problem.verbose);
            fas->atol(1e-8);
            fas->rtol(1e-10);
            fas->max_it(100);

            fas->set_functions(level_functions);
            fas->solve(x);

        }



        void RMTR_test()
        {
            // intial guess
            DVectord x = values(problem.n_dofs[problem.n_levels -1 ], 0.0);

            std::vector<std::shared_ptr<ExtendedFunction<DSMatrixd, DVectord> > >  level_functions(problem.n_levels);


            for(auto l=0; l < problem.n_levels; l++)
            {
                auto fun = std::make_shared<Bratu1D<DSMatrixd, DVectord> >(problem.n_dofs[l]);
                level_functions[l] = fun;

                // making sure that fine level IG is feasible
                if(l+1 == problem.n_levels)
                    fun->apply_bc_to_initial_guess(x);
            }

            auto tr_strategy_coarse = std::make_shared<utopia::SteihaugToint<DSMatrixd, DVectord, HOMEMADE> >();
            tr_strategy_coarse->atol(1e-9);
            // tr_strategy_coarse->verbose(true);
            auto tr_strategy_fine 	= std::make_shared<utopia::SteihaugToint<DSMatrixd, DVectord, HOMEMADE> >();
            tr_strategy_fine->atol(1e-9);
            // tr_strategy_fine->verbose(true);

            // we should apply BC conditions in symmetric way... ASAP...
            tr_strategy_coarse->set_preconditioner(std::make_shared<IdentityPreconditioner<DVectord> >());
            // tr_strategy_fine->set_preconditioner(std::make_shared<InvDiagPreconditioner<DSMatrixd, DVectord> >());

            tr_strategy_coarse->set_preconditioner(std::make_shared<IdentityPreconditioner<DVectord> >());
            // tr_strategy_fine->set_preconditioner(std::make_shared<InvDiagPreconditioner<DSMatrixd, DVectord> >());


            // auto rmtr = std::make_shared<RMTR<DSMatrixd, DVectord, SECOND_ORDER>  >(tr_strategy_coarse, tr_strategy_fine);
            auto rmtr = std::make_shared<RMTR<DSMatrixd, DVectord, FIRST_ORDER>  >(problem.n_levels);
            rmtr->set_transfer_operators(problem.prolongations, problem.restrictions);

            rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
            rmtr->set_fine_tr_strategy(tr_strategy_fine);


            rmtr->max_it(50);
            rmtr->max_coarse_it(1);
            rmtr->max_smoothing_it(3);
            rmtr->delta0(10);
            rmtr->atol(1e-5);
            rmtr->rtol(1e-10);
            rmtr->set_grad_smoothess_termination(0.000001);
            rmtr->set_eps_grad_termination(1e-7);

            rmtr->verbose(problem.verbose);
            rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);
            rmtr->set_functions(level_functions);


            rmtr->solve(x);
        }



        void RMTR_inf_test()
        {
            // intial guess
            DVectord x = values(problem.n_dofs[problem.n_levels -1 ], 0.0);

            std::vector<std::shared_ptr<ExtendedFunction<DSMatrixd, DVectord> > >  level_functions(problem.n_levels);
            for(auto l=0; l < problem.n_levels; l++)
            {
                auto fun = std::make_shared<Bratu1D<DSMatrixd, DVectord> >(problem.n_dofs[l]);
                level_functions[l] = fun;

                // making sure that fine level IG is feasible
                if(l+1 == problem.n_levels)
                    fun->apply_bc_to_initial_guess(x);
            }


            auto lsolver = std::make_shared<LUDecomposition<DSMatrixd, DVectord> >();
            auto tr_strategy_fine =  std::make_shared<utopia::TaoQPSolver<DSMatrixd, DVectord> >(lsolver);
            tr_strategy_fine->set_linear_solver(std::make_shared<GMRES<DSMatrixd, DVectord>>(PCJACOBI));

            auto tr_strategy_coarse =  std::make_shared<utopia::TaoQPSolver<DSMatrixd, DVectord> >(lsolver);
            tr_strategy_coarse->set_linear_solver(std::make_shared<GMRES<DSMatrixd, DVectord>>(PCLU));

            auto rmtr = std::make_shared<RMTR_inf<DSMatrixd, DVectord, SECOND_ORDER>  >(problem.n_levels);
            rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
            rmtr->set_fine_tr_strategy(tr_strategy_fine);


            rmtr->set_transfer_operators(problem.prolongations, problem.restrictions);

            rmtr->max_it(1000);
            rmtr->max_coarse_it(1);
            rmtr->max_smoothing_it(1);
            rmtr->delta0(1);
            rmtr->atol(1e-5);
            rmtr->rtol(1e-10);
            rmtr->set_grad_smoothess_termination(0.000001);
            rmtr->set_eps_grad_termination(1e-7);

            rmtr->verbose(problem.verbose);
            // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
            rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);

            rmtr->set_functions(level_functions);
            rmtr->solve(x);
        }



        void RMTR_inf_bound_test()
        {
            // intial guess
            DVectord x = values(problem.n_dofs[problem.n_levels -1 ], 0.0);

            // upper, lower bound...
            DVectord ub, lb;
            std::vector<std::shared_ptr<ExtendedFunction<DSMatrixd, DVectord> > >  level_functions(problem.n_levels);
            for(auto l=0; l < problem.n_levels; l++)
            {
                auto fun = std::make_shared<Bratu1D<DSMatrixd, DVectord> >(problem.n_dofs[l]);
                level_functions[l] = fun;

                // making sure that fine level IG is feasible
                if(l+1 == problem.n_levels)
                {
                    fun->apply_bc_to_initial_guess(x);
                    fun->generate_constraints(lb, ub, -10, 0.1);
                }
            }


            // Utopia::instance().set("log_output_path", "benchmark.csv");

            auto lsolver = std::make_shared<LUDecomposition<DSMatrixd, DVectord> >();
            auto tr_strategy_fine =  std::make_shared<utopia::TaoQPSolver<DSMatrixd, DVectord> >(lsolver);
            tr_strategy_fine->set_linear_solver(std::make_shared<GMRES<DSMatrixd, DVectord>>(PCJACOBI));
            tr_strategy_fine->verbose(false);

            auto tr_strategy_coarse =  std::make_shared<utopia::TaoQPSolver<DSMatrixd, DVectord> >(lsolver);
            tr_strategy_coarse->set_linear_solver(std::make_shared<GMRES<DSMatrixd, DVectord>>(PCLU));
            // tr_strategy_coarse->verbose(true);
            tr_strategy_coarse->verbose(false);

            auto rmtr = std::make_shared<RMTR_inf<DSMatrixd, DVectord, SECOND_ORDER>  >(problem.n_levels);
            rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
            rmtr->set_fine_tr_strategy(tr_strategy_fine);

            rmtr->set_transfer_operators(problem.prolongations, problem.restrictions);

            rmtr->max_it(1000);
            rmtr->max_coarse_it(1);
            rmtr->max_smoothing_it(1);
            rmtr->delta0(1);
            rmtr->atol(1e-5);
            rmtr->rtol(1e-10);
            rmtr->set_grad_smoothess_termination(0.000001);
            rmtr->set_eps_grad_termination(1e-7);

            rmtr->verbose(problem.verbose);
            // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
            rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);

            rmtr->set_functions(level_functions);

               auto box = make_box_constaints(make_ref(lb), make_ref(ub));
            rmtr->set_box_constraints(box);
            rmtr->solve(x);
        }

    private:

        BratuMultilevelTestProblem<DSMatrixd, DVectord> problem;

    };

#endif //WITH_PETSC


    void runNonlinearMultilevelSolverTest()
    {
        UTOPIA_UNIT_TEST_BEGIN("NonlinearMultilevelSolverTest");
        #ifdef  WITH_PETSC
            NonlinearBratuSolverTest(3, true, false).run();
        #endif
        UTOPIA_UNIT_TEST_END("NonlinearMultilevelSolverTest");

    }
}
