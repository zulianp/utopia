#include "utopia.hpp"
#include "utopia_SolverTest.hpp"
#include "test_problems/utopia_TestProblems.hpp"

namespace utopia
{


#ifdef  WITH_PETSC
	class NonlinearBratuSolverTest {
	public:

		typedef UTOPIA_SIZE_TYPE(DVectord) SizeType;
		typedef UTOPIA_SCALAR(DVectord) Scalar;

		NonlinearBratuSolverTest(const SizeType & n_levels = 2): 
					n_coarse_(10), 
					n_levels_(n_levels) 
		{ 

			assert(n_coarse_ > 0);
			assert(n_levels_ > 1);

			n_dofs_.resize(n_levels_);
			n_dofs_[0] = n_coarse_;

			for(SizeType i = 1; i < n_levels_; ++i) 
				n_dofs_[i] = n_dofs_[i-1] * 2;


			prolongations_.resize(n_levels_ - 1);
			restrictions_.resize(n_levels_ - 1);

			for(SizeType i = 0; i < n_levels_ - 1; ++i) {
				const auto n_coarse = n_dofs_[i];
				const auto n_fine   = n_dofs_[i + 1];
				prolongations_[i] = std::make_shared<DSMatrixd>(sparse(n_fine, n_coarse, 2));
				auto &I = *prolongations_[i];

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

			// remove BC contributions 
			{
				auto &I = *prolongations_.back();

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
			// very not nice solution... 
			for(SizeType i = 0; i < prolongations_.size(); ++i) 
			{
				auto &I = *prolongations_[i];
				DSMatrixd R =  0.25 * transpose(I); 
				restrictions_[i] = std::make_shared<DSMatrixd>(R);
			}
		}		
		
		void run()
		{
			UTOPIA_RUN_TEST(TR_test); 
			UTOPIA_RUN_TEST(TR_constraint_test); 
			UTOPIA_RUN_TEST(newton_MG_test); 
			UTOPIA_RUN_TEST(NMG_test); 
			UTOPIA_RUN_TEST(RMTR_L2_test); 
		}


	    void TR_test()
	    {
	    	Bratu1D<DSMatrixd, DVectord> fun(n_coarse_); 
	    	DVectord x = values(n_coarse_, 1.0); 
	    	fun.apply_bc_to_initial_guess(x); 

	    	Parameters params;
			params.atol(1e-10);
			params.rtol(1e-10);
			params.stol(1e-10);
			params.verbose(true);
			
			auto subproblem = std::make_shared<utopia::KSP_TR<DSMatrixd, DVectord> >();
			TrustRegion<DSMatrixd, DVectord> tr_solver(subproblem);
			tr_solver.set_parameters(params);
			tr_solver.solve(fun, x);
	    }


	    void TR_constraint_test()
	    {
	    	Bratu1D<DSMatrixd, DVectord> fun(n_coarse_); 
	    	DVectord x = values(n_coarse_, 1.0); 
	    	fun.apply_bc_to_initial_guess(x); 

	    	DVectord ub, lb; 
	    	fun.generate_constraints(lb, ub); 
	    	auto box = make_box_constaints(make_ref(lb), make_ref(ub)); 

	    	Parameters params;
			params.atol(1e-6);
			params.rtol(1e-10);
			params.stol(1e-10);
			params.verbose(true);

	        auto lsolver = std::make_shared<LUDecomposition<DSMatrixd, DVectord> >();
	        auto qp_solver = std::make_shared<TaoTRSubproblem<DSMatrixd, DVectord> >(lsolver); 

	        TrustRegionVariableBound<DSMatrixd, DVectord>  tr_solver(qp_solver); 
	        tr_solver.set_box_constraints(box); 
			tr_solver.set_parameters(params);
			tr_solver.solve(fun, x);
	    }	    



	    void newton_MG_test()
	    {
	    	Bratu1D<DSMatrixd, DVectord> fun(n_dofs_[n_levels_ - 1]); 
	    	DVectord x = values(n_dofs_[n_levels_ - 1], 1.0); 
	    	fun.apply_bc_to_initial_guess(x); 

	    	auto lsolver = std::make_shared<utopia::BiCGStab<DSMatrixd, DVectord> >();
            Newton<utopia::DSMatrixd, utopia::DVectord> newton(lsolver);
            
            auto direct_solver = std::make_shared<LUDecomposition<DSMatrixd, DVectord>>();
            auto gs = std::make_shared<GaussSeidel<DSMatrixd, DVectord> >();
            auto multigrid = std::make_shared<Multigrid<DSMatrixd, DVectord>  >(gs, direct_solver);
            
            multigrid->set_transfer_operators(prolongations_);
            multigrid->must_generate_masks(false); 
            multigrid->verbose(false);
            
            newton.set_linear_solver(multigrid);
            newton.verbose(true); 
            newton.atol(1e-9);
            newton.rtol(1e-10);
            newton.solve(fun, x);
	    }	    





	    void NMG_test()
	    {
	    	// this eq. is known to have problems without globalization... 
	    	std::vector<std::shared_ptr<ExtendedFunction<DSMatrixd, DVectord> > >  level_functions(n_levels_); 

	    	for(auto l=0; l < n_levels_-1; l++)
	    	{
		    	Bratu1D<DSMatrixd, DVectord> fun(n_dofs_[l]); 
		    	level_functions[l] = std::make_shared<Bratu1D<DSMatrixd, DVectord> >(fun); 
		    }
	        
	        auto direct_solver = std::make_shared<LUDecomposition<DSMatrixd, DVectord>>();
	        auto coarse_solver = std::make_shared<Newton<utopia::DSMatrixd, utopia::DVectord> >(direct_solver);
	       	auto strategy = std::make_shared<utopia::Backtracking<utopia::DSMatrixd, utopia::DVectord>>();
	        coarse_solver->set_line_search_strategy(strategy);
	        coarse_solver->atol(1e-9); 
	        coarse_solver->max_it(1); 

	        // subject to change 
	        // auto smoother = std::make_shared<NonLinearJacobi<DSMatrixd, DVectord> >();
	        auto smoother = std::make_shared<NonLinearGMRES<DSMatrixd, DVectord> >();

	        auto fas = std::make_shared<NonLinearMultigrid<DSMatrixd, DVectord>  >(smoother, coarse_solver);
	        fas->set_transfer_operators(prolongations_, restrictions_);

			fas->pre_smoothing_steps(5); 
        	fas->post_smoothing_steps(5); 
	        fas->verbose(true); 
	        fas->atol(1e-8); 
	        fas->rtol(1e-10);
	        fas->max_it(50);

	        fas->set_functions(level_functions); 
	        
	        Bratu1D<DSMatrixd, DVectord> fun_fine(n_dofs_[n_levels_-1]); 
	        DVectord x = values(n_dofs_[n_levels_ -1 ], 0.0); 
	        fun_fine.apply_bc_to_initial_guess(x); 

	        fas->solve(fun_fine, x); 
	    }	 



	    void RMTR_L2_test()
	    {
	    	// this eq. is known to have problems without globalization... 
	    	std::vector<std::shared_ptr<ExtendedFunction<DSMatrixd, DVectord> > >  level_functions(n_levels_); 

	    	for(auto l=0; l < n_levels_-1; l++)
	    	{
		    	Bratu1D<DSMatrixd, DVectord> fun(n_dofs_[l]); 
		    	level_functions[l] = std::make_shared<Bratu1D<DSMatrixd, DVectord> >(fun); 
		    }
	        
	        auto tr_strategy_coarse = std::make_shared<utopia::KSP_TR<DSMatrixd, DVectord> >("gltr");
	        tr_strategy_coarse->atol(1e-12); 
	        tr_strategy_coarse->rtol(1e-12); 
	        tr_strategy_coarse->pc_type("lu"); 


	        auto tr_strategy_fine = std::make_shared<utopia::KSP_TR<DSMatrixd, DVectord> >("gltr");
	        tr_strategy_fine->atol(1e-15); 
	        tr_strategy_fine->rtol(1e-15); 
	        tr_strategy_fine->pc_type("jacobi");   

        	auto rmtr = std::make_shared<RMTR<DSMatrixd, DVectord>  >(tr_strategy_coarse, tr_strategy_fine);
	        rmtr->set_transfer_operators(prolongations_, restrictions_);

	        rmtr->max_it(1000); 
	        rmtr->set_max_coarse_it(1); 
	        rmtr->set_max_smoothing_it(1); 
	        rmtr->delta0(1); 
	        rmtr->atol(1e-6); 
	        rmtr->set_grad_smoothess_termination(0.000001); 
			
			rmtr->verbose(true); 
			// rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE); 
			rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL); 

	        rmtr->set_functions(level_functions); 
	        
	        Bratu1D<DSMatrixd, DVectord> fun_fine(n_dofs_[n_levels_-1]); 
	        DVectord x = values(n_dofs_[n_levels_ -1 ], 1.0); 
	        fun_fine.apply_bc_to_initial_guess(x); 

	        rmtr->solve(fun_fine, x); 
	    }	 






		
	private:
		SizeType n_coarse_;
		SizeType n_levels_; 
		std::vector<SizeType> n_dofs_;

		std::vector<std::shared_ptr<DSMatrixd>> prolongations_;
		std::vector<std::shared_ptr<DSMatrixd>> restrictions_;

	};

#endif //WITH_PETSC


	void runNonlinearMultilevelSolverTest()
	{

		UTOPIA_UNIT_TEST_BEGIN("runNonlinearMultilevelSolverTest");
		#ifdef  WITH_PETSC
			NonlinearBratuSolverTest(3).run();
		#endif		
		UTOPIA_UNIT_TEST_END("runNonlinearMultilevelSolverTest");				
	}
}
