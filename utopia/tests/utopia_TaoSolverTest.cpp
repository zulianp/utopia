#include "utopia_TaoSolverTest.hpp"
#include "utopia.hpp"

#include "test_problems/utopia_TestProblems.hpp"
#include "utopia_petsc_TaoSolver.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"
#include "utopia_QuadraticFunction.hpp"

#ifdef WITH_PETSC

namespace utopia {
	void petsc_tao_solve_simple()
	{	
		TestFunctionND_1<DMatrixd, DVectord> fun(10);
		TaoSolver<DMatrixd, DVectord> tao(std::make_shared<Factorization<DMatrixd, DVectord>>());
		DVectord x = zeros(10);
		tao.solve(fun, x);

		DVectord expected = values(10, 0.468919);
		assert(approxeq(x, expected));
	}

	void petsc_tao_solve_vi()
	{
		const SizeType n = 10;

		DSMatrixd m = zeros(n, n);
		assemble_laplacian_1D(n, m);

		{
			Range r = row_range(m);
			Write<DSMatrixd> w(m);
			if(r.begin() == 0) {
				m.set(0, 0, 1.);
				m.set(0, 1, 0);
			}

			if(r.end() == n) {
				m.set(n-1, n-1, 1.);
				m.set(n-1, n-2, 0);
			}
		}

		DVectord rhs = values(n, 1.);
		{ 
			Range r = range(rhs);
			Write<DVectord> w(rhs);

			if(r.inside(0)) {
				rhs.set(0, 0);
			}

			if(r.inside(n-1)) {
				rhs.set(n-1, 0.);
			}
		}

		DVectord upper_bound = values(n, 100.0);
		DVectord x = zeros(n);

		auto box = make_upper_bound_constraints(make_ref(upper_bound));

		QuadraticFunction<DSMatrixd, DVectord> fun(make_ref(m), make_ref(rhs));
		TaoSolver<DSMatrixd, DVectord> tao(std::make_shared<Factorization<DSMatrixd, DVectord>>());
		// tao.set_box_constraints(box);
		tao.solve(fun, x);
	}

	void run_tao_solver_test()
	{
		UTOPIA_RUN_TEST(petsc_tao_solve_simple);
		// UTOPIA_RUN_TEST(petsc_tao_solve_vi);
	}
}

#else
namespace utopia { void run_tao_solver_test() {} }
#endif //WITH_PETSC
