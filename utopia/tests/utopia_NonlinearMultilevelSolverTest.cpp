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
			assert(problem.n_coarse > 0);
			assert(problem.n_levels > 1);

			problem.n_dofs.resize(problem.n_levels);
			problem.n_dofs[0] = problem.n_coarse;

			for(SizeType i = 1; i < problem.n_levels; ++i)
				problem.n_dofs[i] = (problem.n_dofs[i-1] - 1) * 2 + 1;


			problem.prolongations.resize(problem.n_levels - 1);
			problem.restrictions.resize(problem.n_levels - 1);

			for(SizeType i = 0; i < problem.n_levels - 1; ++i) {
				const auto n_coarse = problem.n_dofs[i];
				const auto n_fine   = problem.n_dofs[i + 1];
				problem.prolongations[i] = std::make_shared<DSMatrixd>(sparse(n_fine, n_coarse, 2));
				auto &I = *problem.prolongations[i];

				Write<DSMatrixd> w_(I);
				auto r = row_range(I);

				SizeType j = r.begin()/2;
				for(auto k = r.begin(); k < r.end(); k += 2, ++j) {
					I.set(k, j, 1.);

					if(j + 1 < n_coarse) {
						I.set(k + 1, j, 0.5);
						I.set(k + 1, j + 1, 0.5);
					}
				}
			}

			if(problem.remove_BC_contributions)
			{
				auto &I = *problem.prolongations.back();

				Write<DSMatrixd> w_(I);
				auto rr = row_range(I);

				if(rr.inside(0)) {
					I.set(0, 0, 0.);
				}

				auto last_node_h = size(I).get(0) - 1;
				auto last_node_H = size(I).get(1) - 1;
				if(rr.inside(last_node_h)) {
					I.set(last_node_h, last_node_H, 0.);
				}
			}

			// restrictions, but let's use them as projections...
			// not very nice solution, but I am lazy to do something more sophisticated just for testing purposes...
			for(SizeType i = 0; i < problem.prolongations.size(); ++i)
			{
				auto &I = *problem.prolongations[i];
				DSMatrixd R =  0.5*  transpose(I);
				problem.restrictions[i] = std::make_shared<DSMatrixd>(R);
			}
		}

		void run()
		{
			UTOPIA_RUN_TEST(TR_test);
			UTOPIA_RUN_TEST(TR_constraint_test);

			UTOPIA_RUN_TEST(newton_MG_test);
			UTOPIA_RUN_TEST(FAS_test);

			UTOPIA_RUN_TEST(RMTR_test);
			UTOPIA_RUN_TEST(RMTR_inf_test);

			UTOPIA_RUN_TEST(RMTR_inf_bound_test);
		}


	    void TR_test()
	    {
	    	Bratu1D<DSMatrixd, DVectord> fun(problem.n_coarse);
	    	DVectord x = values(problem.n_coarse, 1.0);
	    	fun.apply_bc_to_initial_guess(x);

	    	Parameters params;
			params.atol(1e-10);
			params.rtol(1e-10);
			params.stol(1e-10);
			params.verbose(problem.verbose);

			auto subproblem = std::make_shared<utopia::KSP_TR<DSMatrixd, DVectord> >();
			TrustRegion<DSMatrixd, DVectord> tr_solver(subproblem);
			tr_solver.set_parameters(params);
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

	    	Parameters params;
			params.atol(1e-6);
			params.rtol(1e-10);
			params.stol(1e-10);
			params.verbose(problem.verbose);

	        auto lsolver = std::make_shared<LUDecomposition<DSMatrixd, DVectord> >();
	        auto qp_solver = std::make_shared<TaoTRSubproblem<DSMatrixd, DVectord> >(lsolver);

	        TrustRegionVariableBound<DSMatrixd, DVectord>  tr_solver(qp_solver);
	        tr_solver.set_box_constraints(box);
			tr_solver.set_parameters(params);
			tr_solver.solve(fun, x);
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
            multigrid->verbose(problem.verbose);
            multigrid->atol(1e-11);

            newton.set_linear_solver(multigrid);
            newton.verbose(problem.verbose);
            newton.atol(1e-9);
            newton.rtol(1e-10);
            newton.solve(fun, x);
	    }


	    void FAS_test()
	    {
	    	if(mpi_world_size() > 1)
	    		return;

	    	// intial guess
	        DVectord x = values(problem.n_dofs[problem.n_levels -1 ], 0.0);

	    	std::vector<std::shared_ptr<ExtendedFunction<DSMatrixd, DVectord> > >  level_functions(problem.n_levels);


	    	for(auto l=0; l < problem.n_levels; l++)
	    	{
		    	Bratu1D<DSMatrixd, DVectord> fun(problem.n_dofs[l]);
		    	level_functions[l] = std::make_shared<Bratu1D<DSMatrixd, DVectord> >(fun);

		    	// making sure that fine level IG is feasible
		    	if(l+1 == problem.n_levels)
		    		fun.apply_bc_to_initial_guess(x);
		    }


	        auto direct_solver = std::make_shared<LUDecomposition<DSMatrixd, DVectord>>();
	        auto coarse_solver = std::make_shared<Newton<utopia::DSMatrixd, utopia::DVectord> >(direct_solver);
	       	auto strategy = std::make_shared<utopia::Backtracking<utopia::DSMatrixd, utopia::DVectord>>();
	        coarse_solver->set_line_search_strategy(strategy);
	        coarse_solver->atol(1e-7);
	        coarse_solver->max_it(1);

	        // subject to change
	        auto smoother = std::make_shared<NonLinearJacobi<DSMatrixd, DVectord> >();
	        smoother->damping_parameter(0.3);
	        // auto smoother = std::make_shared<NonLinearGMRES<DSMatrixd, DVectord> >();


	        auto fas = std::make_shared<FAS<DSMatrixd, DVectord>  >(smoother, coarse_solver);
	        fas->set_transfer_operators(problem.prolongations, problem.restrictions, problem.restrictions);

			fas->pre_smoothing_steps(3);
        	fas->post_smoothing_steps(3);
	        fas->verbose(problem.verbose);
	        fas->atol(1e-8);
	        fas->rtol(1e-10);
	        fas->max_it(10);

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
		    	Bratu1D<DSMatrixd, DVectord> fun(problem.n_dofs[l]);
		    	level_functions[l] = std::make_shared<Bratu1D<DSMatrixd, DVectord> >(fun);

		    	// making sure that fine level IG is feasible
		    	if(l+1 == problem.n_levels)
		    		fun.apply_bc_to_initial_guess(x);
		    }

	        auto tr_strategy_coarse = std::make_shared<utopia::SteihaugToint<DSMatrixd, DVectord, HOMEMADE> >();
	        auto tr_strategy_fine 	= std::make_shared<utopia::SteihaugToint<DSMatrixd, DVectord, HOMEMADE> >();

        	// auto rmtr = std::make_shared<RMTR<DSMatrixd, DVectord, SECOND_ORDER>  >(tr_strategy_coarse, tr_strategy_fine);
        	auto rmtr = std::make_shared<RMTR<DSMatrixd, DVectord, GALERKIN>  >(tr_strategy_coarse, tr_strategy_fine);
	        rmtr->set_transfer_operators(problem.prolongations, problem.restrictions);

	        rmtr->max_it(1000);
	        rmtr->max_coarse_it(1);
	        rmtr->max_smoothing_it(1);
	        rmtr->delta0(1);
	        rmtr->atol(1e-6);
	        rmtr->rtol(1e-10);
	        rmtr->set_grad_smoothess_termination(0.000001);
	        rmtr->set_eps_grad_termination(1e-7);

			rmtr->verbose(problem.verbose);
			// rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
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
		    	Bratu1D<DSMatrixd, DVectord> fun(problem.n_dofs[l]);
		    	level_functions[l] = std::make_shared<Bratu1D<DSMatrixd, DVectord> >(fun);

		    	// making sure that fine level IG is feasible
		    	if(l+1 == problem.n_levels)
		    		fun.apply_bc_to_initial_guess(x);
		    }


		    auto lsolver = std::make_shared<LUDecomposition<DSMatrixd, DVectord> >();
        	auto tr_strategy_fine = std::make_shared<TaoTRSubproblem<DSMatrixd, DVectord> >(lsolver);
        	tr_strategy_fine->pc_type("jacobi");

        	auto tr_strategy_coarse = std::make_shared<TaoTRSubproblem<DSMatrixd, DVectord> >(lsolver);
        	tr_strategy_coarse->pc_type("lu");

        	auto rmtr = std::make_shared<RMTR_inf<DSMatrixd, DVectord, SECOND_ORDER>  >(tr_strategy_coarse, tr_strategy_fine);
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
		    	Bratu1D<DSMatrixd, DVectord> fun(problem.n_dofs[l]);
		    	level_functions[l] = std::make_shared<Bratu1D<DSMatrixd, DVectord> >(fun);

		    	// making sure that fine level IG is feasible
		    	if(l+1 == problem.n_levels)
		    	{
		    		fun.apply_bc_to_initial_guess(x);
		    		fun.generate_constraints(lb, ub, -10, 0.1);
		    	}
		    }


		    // Utopia::instance().set("log_output_path", "benchmark.csv");

		    auto lsolver = std::make_shared<LUDecomposition<DSMatrixd, DVectord> >();
        	auto tr_strategy_fine = std::make_shared<TaoTRSubproblem<DSMatrixd, DVectord> >(lsolver);
        	tr_strategy_fine->pc_type("jacobi");
        	tr_strategy_fine->verbose(problem.verbose);

        	auto tr_strategy_coarse = std::make_shared<TaoTRSubproblem<DSMatrixd, DVectord> >(lsolver);
        	tr_strategy_coarse->pc_type("lu");
        	// tr_strategy_coarse->verbose(true);
        	tr_strategy_coarse->verbose(problem.verbose);

        	auto rmtr = std::make_shared<RMTR_inf<DSMatrixd, DVectord, SECOND_ORDER>  >(tr_strategy_coarse, tr_strategy_fine);
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
		UTOPIA_UNIT_TEST_BEGIN("runNonlinearMultilevelSolverTest");
		#ifdef  WITH_PETSC
			NonlinearBratuSolverTest(4, true, false).run();
		#endif
		UTOPIA_UNIT_TEST_END("runNonlinearMultilevelSolverTest");

	}
}
