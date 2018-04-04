#include "utopia_TaoSolverTest.hpp"
#include "utopia.hpp"

#include "test_problems/utopia_TestProblems.hpp"
#include "utopia_petsc_TaoSolver.hpp"

#ifdef WITH_PETSC

namespace utopia {
	void run_simple_fun()
	{	
		TestFunctionND_1<DMatrixd, DVectord> fun(10);
		TaoSolver<DMatrixd, DVectord> tao(std::make_shared<Factorization<DMatrixd, DVectord>>());
		DVectord x = zeros(10);
		tao.solve(fun, x);

		DVectord expected = values(10, 0.468919);
		assert(approxeq(x, expected));
	}

	void run_tao_solver_test()
	{
		UTOPIA_RUN_TEST(run_simple_fun);
	}
}

#else
namespace utopia { void run_tao_solver_test() {} }
#endif //WITH_PETSC
