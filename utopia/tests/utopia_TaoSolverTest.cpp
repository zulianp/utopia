#include "utopia_TaoSolverTest.hpp"
#include "utopia.hpp"

#ifdef WITH_PETSC

#include "test_problems/utopia_TestProblems.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"
#include "utopia_QuadraticFunction.hpp"
#include "utopia_petsc_TaoSolver.hpp"

namespace utopia {
	void petsc_tao_solve_simple()
	{	
		TestFunctionND_1<DMatrixd, DVectord> fun(10);
		TaoSolver<DMatrixd, DVectord> tao(std::make_shared<Factorization<DMatrixd, DVectord>>());
		DVectord x = zeros(10);
		tao.set_type("blmvm");
		tao.solve(fun, x);

		DVectord expected = values(10, 0.468919);
		utopia_test_assert(approxeq(x, expected));
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
		utopia_test_assert(n_diff < 1e-10);


	}

	void petsc_tao_solve_mg()
	{
		DVectord rhs;
		DSMatrixd A, I_1, I_2, I_3;
		
		const std::string data_path = Utopia::instance().get("data_path");
		
		read(data_path + "/laplace/matrices_for_petsc/f_rhs", rhs);
		read(data_path + "/laplace/matrices_for_petsc/f_A", A);
		read(data_path + "/laplace/matrices_for_petsc/I_2", I_2);
		read(data_path + "/laplace/matrices_for_petsc/I_3", I_3);
		
		std::vector<std::shared_ptr<DSMatrixd>> interpolation_operators;
		interpolation_operators.push_back(make_ref(I_2));
		interpolation_operators.push_back(make_ref(I_3));
		
		auto smoother      = std::make_shared<GaussSeidel<DSMatrixd, DVectord>>();
		auto linear_solver = std::make_shared<ConjugateGradient<DSMatrixd, DVectord>>();
		Multigrid<DSMatrixd, DVectord> multigrid(smoother, linear_solver);
		multigrid.set_transfer_operators(std::move(interpolation_operators));
		DVectord x = zeros(A.size().get(0));
		// DVectord upper_bound = values(A.size().get(0), 0.003);
		// auto box = make_upper_bound_constraints(make_ref(upper_bound));

		QuadraticFunction<DSMatrixd, DVectord> fun(make_ref(A), make_ref(rhs));
		TaoSolver<DSMatrixd, DVectord> tao(make_ref(multigrid));
		
		// multigrid.verbose(true);
		multigrid.max_it(20);
		multigrid.atol(1e-15);
		multigrid.stol(1e-15);
		multigrid.rtol(1e-15);
		//constraints do not work with mg because system contains lagr mult
		// tao.set_box_constraints(box);
		tao.solve(fun, x);
	}

	void petsc_tao_tr_bound()
	{
		const SizeType n = 100;

		DSMatrixd m;
		DVectord rhs, upper_bound;
		ExampleTestCase2<DSMatrixd, DVectord> example;
		example.getOperators(n, m, rhs, upper_bound);
		DVectord x = zeros(n);

		const double scale_factor = 1e-10;
		rhs *= scale_factor;
		upper_bound *= scale_factor;

		auto box = make_upper_bound_constraints(make_ref(upper_bound));
		QuadraticFunction<DSMatrixd, DVectord> fun(make_ref(m), make_ref(rhs));

		// auto lsolver = std::make_shared<LUDecomposition<DSMatrixd, DVectord> >();
		auto lsolver = std::make_shared<BiCGStab<DSMatrixd, DVectord> >();
        auto qp_solver = std::make_shared<TaoTRSubproblem<DSMatrixd, DVectord> >(lsolver); 

        TrustRegionVariableBound<DSMatrixd, DVectord>  tr_solver(qp_solver); 
        tr_solver.set_box_constraints(box); 
        tr_solver.verbose(false); 
        tr_solver.solve(fun, x); 

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
		utopia_test_assert(n_diff < 1e-10);

	}

	void run_tao_solver_test()
	{
		UTOPIA_UNIT_TEST_BEGIN("PetscTaoTest");
		//does not work yet missing ksp for dense matrix
		// UTOPIA_RUN_TEST(petsc_tao_solve_simple);
		// UTOPIA_RUN_TEST(petsc_tao_solve_vi);
		UTOPIA_RUN_TEST(petsc_tao_solve_mg);
		UTOPIA_RUN_TEST(petsc_tao_tr_bound);
		UTOPIA_UNIT_TEST_END("PetscTaoTest");
	}
}

#else
namespace utopia { void run_tao_solver_test() {} }
#endif //WITH_PETSC
