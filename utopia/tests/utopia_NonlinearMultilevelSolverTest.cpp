#include "utopia.hpp"
#include "utopia_SolverTest.hpp"
#include "test_problems/utopia_TestProblems.hpp"

namespace utopia
{


#ifdef  WITH_PETSC
	class NonlinearMultilevelSolverTest {
	public:
		
		void run()
		{
			UTOPIA_RUN_TEST(TR_Bratu_test); 
		}


	    void TR_Bratu_test()
	    {
	    	Bratu1D<DSMatrixd, DVectord> fun(_n); 
	    	DVectord x = values(_n, 1.0); 
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


		NonlinearMultilevelSolverTest()
		: _n(10) { }
		
	private:
		int _n;
	};

#endif //WITH_PETSC


	void runNonlinearMultilevelSolverTest()
	{

		UTOPIA_UNIT_TEST_BEGIN("runNonlinearMultilevelSolverTest");
		#ifdef  WITH_PETSC
				NonlinearMultilevelSolverTest().run();
		#endif		
		UTOPIA_UNIT_TEST_END("runNonlinearMultilevelSolverTest");				
	}
}
