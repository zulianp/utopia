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
		const SizeType n = 50;

		DSMatrixd m;
		DVectord rhs, upper_bound;
		ExampleTestCase2<DSMatrixd, DVectord> example;
		example.getOperators(n, m, rhs, upper_bound);
		DVectord x = zeros(n);

		//-tao_type blmvm, gpcg seem to deliver the correct result
		//TODO try lcl with TaoSetStateIS
		//-tao_type test is usefull for checking correctness of gradient
		const double scale_factor = 1e-8;
		m *= scale_factor;
		rhs *= scale_factor;
		upper_bound *= scale_factor;

		auto box = make_upper_bound_constraints(make_ref(upper_bound));

		QuadraticFunction<DSMatrixd, DVectord> fun(make_ref(m), make_ref(rhs));
		TaoSolver<DSMatrixd, DVectord> tao(std::make_shared<Factorization<DSMatrixd, DVectord>>());
		tao.set_box_constraints(box);
		tao.solve(fun, x);

		// disp(x);
		x *= 1./scale_factor;
		// x.implementation().set_name("v");
		// write("x.m", x);
	}

	void run_tao_solver_test()
	{
		UTOPIA_RUN_TEST(petsc_tao_solve_simple);
		UTOPIA_RUN_TEST(petsc_tao_solve_vi);
	}
}

#else
namespace utopia { void run_tao_solver_test() {} }
#endif //WITH_PETSC
