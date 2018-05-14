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
					n_coarse_(5), 
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

			// restrictions
			for(SizeType i = 0; i < prolongations_.size(); ++i) 
			{
				auto &I = *prolongations_[i];
				DSMatrixd R = 1.0/2.0 * transpose(I); 
				restrictions_[i] = std::make_shared<DSMatrixd>(R);
			}





		}		
		
		void run()
		{
			UTOPIA_RUN_TEST(TR_test); 
			UTOPIA_RUN_TEST(TR_constraint_test); 
			UTOPIA_RUN_TEST(newton_MG_test); 
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
			params.verbose(false);
			
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
			params.verbose(false);

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
            multigrid->verbose(true);
            
            newton.set_linear_solver(multigrid);
            newton.verbose(true); 
            newton.atol(1e-9);
            newton.rtol(1e-9);
            newton.solve(fun, x);
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
			NonlinearBratuSolverTest().run();
		#endif		
		UTOPIA_UNIT_TEST_END("runNonlinearMultilevelSolverTest");				
	}
}
