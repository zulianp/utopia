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
