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
		tao.set_type("blmvm");
		tao.solve(fun, x);

		DVectord expected = values(10, 0.468919);
		assert(approxeq(x, expected));
	}

	void petsc_tao_solve_vi()
	{
		const SizeType n = 100;

		DSMatrixd m;
		DVectord rhs, upper_bound;
		ExampleTestCase2<DSMatrixd, DVectord> example;
		example.getOperators(n, m, rhs, upper_bound);
		DVectord x = zeros(n);
		auto lsolver = std::make_shared<ConjugateGradient<DSMatrixd, DVectord>>();

		const double scale_factor = 1e-10;
		rhs *= scale_factor;
		upper_bound *= scale_factor;

		auto box = make_upper_bound_constraints(make_ref(upper_bound));

		QuadraticFunction<DSMatrixd, DVectord> fun(make_ref(m), make_ref(rhs));
		TaoSolver<DSMatrixd, DVectord> tao(lsolver);
		// tao.set_ksp_types("bcgs", "jacobi", " ");
		tao.set_box_constraints(box);
		// tao.set_type("tron");
		// tao.set_type("gpcg");
		tao.solve(fun, x);

		x *= 1./scale_factor;

		DVectord xssn = zeros(n);
		SemismoothNewton<DSMatrixd, DVectord, HOMEMADE> ssnewton(std::make_shared<Factorization<DSMatrixd, DVectord>>());
		ssnewton.set_box_constraints(box);
		ssnewton.stol(1e-18);
		ssnewton.atol(1e-18);
		ssnewton.rtol(1e-18);
		ssnewton.solve(m, rhs, xssn);

		xssn *= 1./scale_factor;

		double n_diff = norm2(xssn - x);
		assert(n_diff < 1e-10);


	}

	void run_tao_solver_test()
	{
		UTOPIA_UNIT_TEST_BEGIN("PetscTaoTest");
		//does not work yet missing ksp for dense matrix
		// UTOPIA_RUN_TEST(petsc_tao_solve_simple);
		UTOPIA_RUN_TEST(petsc_tao_solve_vi);
		UTOPIA_UNIT_TEST_END("PetscTaoTest");
	}
}

#else
namespace utopia { void run_tao_solver_test() {} }
#endif //WITH_PETSC
