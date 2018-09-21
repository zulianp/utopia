#include "utopia.hpp"
#include "utopia_SolverTest.hpp"
#include "test_problems/utopia_TestProblems.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"

#include "utopia_ProjectedConjugateGradient.hpp"
#include "utopia_ProjectedGradient.hpp"
#include "utopia_MultiLevelTestProblem.hpp"

namespace utopia {
	/**
	 * @brief      Class to test our nonlinear solvers.
	 *
	 * @todo       add more tests to TR, LS, ...
	 * @todo       different strategies, different LS solve, test all params, nonlinear problems ...
	 *
	 * @tparam     Matrix
	 * @tparam     Vector
	 * @tparam     Scalar
	 */
	template<class Matrix, class Vector, class Scalar>
	class SolverTest {
	public:
		static void print_backend_info()
		{
			if(Utopia::instance().verbose() && mpi_world_rank() == 0) {
				std::cout << "\nBackend: " << backend_info(Vector()).get_name() << std::endl;
			}
		}

		void run()
		{
			print_backend_info();
			UTOPIA_RUN_TEST(newton_cg_test);
			UTOPIA_RUN_TEST(solver_from_params_test);
			UTOPIA_RUN_TEST(tr_test);
			UTOPIA_RUN_TEST(ls_test);
			UTOPIA_RUN_TEST(nl_solve_test);
			UTOPIA_RUN_TEST(dogleg_test);
			UTOPIA_RUN_TEST(st_cg_test); 
			UTOPIA_RUN_TEST(precond_st_cg_test); 
			UTOPIA_RUN_TEST(Quasi_TR_test); 

		}

		class EmptyLSFun : public LeastSquaresFunction<Matrix, Vector> {
		public:
			bool residual(const Vector &/*point*/, Vector &/*result*/) const { return true; }
			bool jacobian(const Vector &/*x*/, Matrix &/*hessian*/) const { return true; }
			bool value(const Vector &, Scalar &val) const { return true; }
			bool update(const Vector &) { return true; }
		};



		void ls_normal_eq()
		{
			LeastSquaresNewton<Matrix, Vector> newton(std::make_shared<ConjugateGradient<Matrix, Vector>>());
			auto ls_strat  = std::make_shared<utopia::Backtracking<Matrix, Vector> >();
			newton.set_line_search_strategy(ls_strat);

			EmptyLSFun fun;
			Vector x0;
			newton.solve(fun, x0);
		}

		void st_cg_test()
        {
            SteihaugToint<Matrix, Vector, HOMEMADE> cg;
            cg.rtol(1e-7);
            cg.atol(1e-6);
            cg.max_it(_n);
            cg.verbose(false);

            Matrix A = sparse(_n, _n, 3);
			assemble_symmetric_laplacian_1D(A, true);

			Vector rhs = values(_n, 975.9);

            {
            	Write<Vector> w(rhs);
            	rhs.set(0, 0.0); 
            	rhs.set(_n-1, 0.0); 
            }			

            Vector x = zeros(size(rhs));

            cg.tr_constrained_solve(A, -1.0 * rhs, x, 1e15);
            utopia_test_assert(approxeq(rhs, A * x, 1e-5));

        }


		void precond_st_cg_test()
        {
            SteihaugToint<Matrix, Vector, HOMEMADE> cg;
            cg.rtol(1e-7);
            cg.atol(1e-6);
            cg.max_it(_n);
            cg.verbose(false);
            cg.set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector> >());
            // cg.set_preconditioner(std::make_shared<IdentityPreconditioner<Matrix, Vector> >());
            

            Matrix A = sparse(_n, _n, 3);
			assemble_symmetric_laplacian_1D(A, true);

			Vector rhs = values(_n, 975.9);

            {
            	Write<Vector> w(rhs);
            	rhs.set(0, 0.0); 
            	rhs.set(_n-1, 0.0); 
            }			

            Vector x = zeros(size(rhs));

            cg.tr_constrained_solve(A, -1.0 * rhs, x, 1e15);
            utopia_test_assert(approxeq(rhs, A * x, 1e-5));
        }



		void nl_solve_test()
		{
			//! [NL solve example]

			//set-up problem
			int n = 10;
			Vector actual = values(n, 2.);
			TestFunctionND_1<Matrix, Vector> fun(n);

			//solve problem
			solve(fun, actual);

			//test outcome...
			Vector expected = values(n, 0.468919);
			utopia_test_assert(approxeq(expected, actual));
			//! [NL solve example]
		}

		void newton_cg_test()
		{
			//! [Newton CG example]
			using namespace std;

			//CG with diagonal preconditioner
			auto linear_solver  = make_shared< ConjugateGradient<Matrix, Vector> >();
			auto preconditioner = make_shared< InvDiagPreconditioner<Matrix, Vector> >();
			linear_solver->set_preconditioner(preconditioner);

			//Newton solver with cg linear solver
			Newton<Matrix, Vector> newton_solver(linear_solver);

			const int n = 10;
			Vector actual   = values(n, 2.);
			Vector expected = values(n, 0.468919);

			TestFunctionND_1<Matrix, Vector> fun(n);

			newton_solver.solve(fun, actual);

			//Check if the result is what we expected
			utopia_test_assert(approxeq(expected, actual));
			//! [Newton CG example]
		}

		void solver_from_params_test()
		{
			Vector x = values(10, 2.0);
			Vector expected = values(x.size().get(0), 0.468919);
			TestFunctionND_1<Matrix, Vector> fun2(x.size().get(0));

			Parameters params;
			params.tol(1e-7);
			params.solver_type(TRUST_REGION_TAG);
			params.lin_solver_type(BICGSTAB_TAG);
			params.linear_solver_verbose(false);
			params.verbose(false);

			// solve(fun2, x, params);

			// monitor test
			// Matrix H = identity(2,2);
			// int i = 0;
			// monitor(i, H);
			// Matrixd blas_mat  = identity(3,3);
			// monitor(i, blas_mat);
		}

		void tr_test()
		{
			// rosenbrock test
			if(mpi_world_size() == 1)
			{
				Vector x = values(10, 2);
				TestFunctionND_1<Matrix, Vector> fun2(x.size().get(0));
				Vector expected = values(x.size().get(0), 0.468919);

				Vector x_w1  = values(4, 10);
				Vector expected_woods = values(4, 1);
				Woods<Matrix, Vector> fun_woods;
				{
					Write<Vector> w1(x_w1);
					x_w1.set(0, -3);
					x_w1.set(1, -1);
					x_w1.set(2, -3);
					x_w1.set(3, -1);
				}


				Parameters params;
				params.atol(1e-10);
				params.rtol(1e-10);
				params.stol(1e-10);
				params.solver_type(TRUST_REGION_TAG);
				params.lin_solver_type(BICGSTAB_TAG);
				params.trust_region_alg(DOGLEG_TAG);
				params.verbose(false);

				// trust_region_solve(fun2, x, params);
				// trust_region_solve(fun_woods, x_w1, params);
				// utopia_test_assert(approxeq(expected, x));


				x = values(10, 2);
				params.trust_region_alg(STEIHAUG_TOINT_TAG);
				trust_region_solve(fun2, x, params);
				utopia_test_assert(approxeq(expected, x));


				// x = values(10, 2);
				// params.trust_region_alg(STEIHAUG_TOINT_TAG);

				// auto strategy_tr = std::make_shared<utopia::SteihaugToint<Matrix, Vector> >();
				// auto precond = std::make_shared< InvDiagPreconditioner<Matrix, Vector> >();
				// strategy_tr->set_preconditioner(precond);
				// strategy_tr->verbose(true);
				// TrustRegion<Matrix, Vector> tr_solver(strategy_tr);

				// tr_solver.verbose(true);
				// tr_solver.solve(fun2, x);
				// utopia_test_assert(approxeq(expected, x));



				x = values(10, 2);
				params.trust_region_alg(CAUCHYPOINT_TAG);
				trust_region_solve(fun2, x, params);
				utopia_test_assert(approxeq(expected, x));

				x = values(10, 2);
				trust_region_solve(fun2, x, params);
				utopia_test_assert(approxeq(expected, x));


				Vector expected_rosenbrock = values(2, 1);

				Rosenbrock<Matrix, Vector> rosenbrock;
				Vector x0 = values(2, 2.0);

				x0 = values(2, 2.0);
				params.trust_region_alg(STEIHAUG_TOINT_TAG);
				params.verbose(false);
				params.atol(1e-13);
				params.rtol(1e-17);
				trust_region_solve(rosenbrock, x0, params);

				auto diff_norm = norm2(expected_rosenbrock - x0);

	            if(diff_norm > 1e-12) {
	                utopia_error("tr_test: STEIHAUG_TOINT_TAG with rosenbrock is failing");
	            }
			}
		}

		void ls_test()
		{
			if(mpi_world_size() == 1) {
				Vector x1 = values(10, 2);
				Vector x2 = values(10, 2);
				TestFunctionND_1<Matrix, Vector> fun2(x1.size().get(0));

				Vector expected = values(x1.size().get(0), 0.468919);

				Parameters params ;
				params.atol(1e-11);
				params.rtol(1e-11);
				params.stol(1e-11);
				params.verbose(false);
				params.linear_solver_verbose(false);
				params.line_search_inner_verbose(false);

				auto lsolver = std::make_shared< ConjugateGradient<Matrix, Vector, HOMEMADE> >();
				Newton<Matrix, Vector> nlsolver1(lsolver);
				Newton<Matrix, Vector> nlsolver2(lsolver);


				auto strategy_sbc = std::make_shared<utopia::SimpleBacktracking<Matrix, Vector> >();
				auto strategy_bc  = std::make_shared<utopia::Backtracking<Matrix, Vector> >();

				strategy_sbc->set_parameters(params);
				strategy_bc->set_parameters(params);

				nlsolver1.set_line_search_strategy(strategy_sbc);
				nlsolver2.set_line_search_strategy(strategy_bc);

				nlsolver1.set_parameters(params);
				nlsolver2.set_parameters(params);


				nlsolver1.solve(fun2, x1);
				nlsolver2.solve(fun2, x2);

				// Woods function test
				Vector x_w1  = values(4, 10);
				Vector x_w2  = values(4, 10);
				Vector expected_woods = values(4, 1);
				{
					Write<Vector> w1(x_w1);
					Write<Vector> w2(x_w2);
					x_w1.set(0, -3);    x_w2.set(0, -3);
					x_w1.set(1, -1);    x_w2.set(1, -1);
					x_w1.set(2, -3);    x_w2.set(2, -3);
					x_w1.set(3, -1);    x_w2.set(3, -1);
				}

				Woods<Matrix, Vector> fun_woods;
				nlsolver1.solve(fun_woods, x_w1);
				nlsolver2.solve(fun_woods, x_w2);

				utopia_test_assert(approxeq(expected_woods, x_w1));
				utopia_test_assert(approxeq(expected_woods, x_w2));

				// rastrigin function test - convergence to local minimum
				Rastrigin<Matrix, Vector> fun_rastrigin;
				Vector x_r1 = values(2, 1), x_r2 = values(2, 1), expected_rastrigin = values(2, 1);
				{
					Write<Vector> w1(x_r1);
					Write<Vector> w2(x_r2);
					Write<Vector> w3(expected_rastrigin);
					x_r1.set(0, -5.12); x_r2.set(0, -5.12);  expected_rastrigin.set(0, -4.97469);
					x_r1.set(1, 5.12); x_r2.set(1, 5.12);  expected_rastrigin.set(1, 4.97469);
				}

				nlsolver1.solve(fun_rastrigin, x_r2);
				nlsolver2.solve(fun_rastrigin, x_r2);

				// rosenbrock test

				Vector expected_rosenbrock = values(2, 1);
				Rosenbrock<Matrix, Vector> rosenbrock_fun;

				Vector x01 = values(2, 2.0), x02 = values(2, 2.0);
				nlsolver1.solve(rosenbrock_fun, x01);
				nlsolver2.solve(rosenbrock_fun, x02);
				utopia_test_assert(approxeq(expected_rosenbrock, x01));
				utopia_test_assert(approxeq(expected_rosenbrock, x02));

				x01 = values(2, 2.0);
				params.verbose(false);
				params.linear_solver_verbose(false);
				line_search_solve(rosenbrock_fun, x01, params);
			}
		}


		void dogleg_test()
		{
			// rosenbrock test
			if(mpi_world_size() == 1)
			{
				Rosenbrock<Matrix, Vector> rosenbrock;
				Vector expected_rosenbrock = values(2, 1);

				auto dogleg = std::make_shared<Dogleg<Matrix, Vector> >();

				Vector x0 = values(2, 2.0);

				TrustRegion<Matrix, Vector> tr_solver(dogleg);
				tr_solver.verbose(false);
				auto cg = std::make_shared<ConjugateGradient<Matrix, Vector> >();
				tr_solver.set_linear_solver(cg);				
				tr_solver.solve(rosenbrock, x0);

				utopia_test_assert(approxeq(expected_rosenbrock, x0));
			}
		}


		void Quasi_TR_test()
		{
			// rosenbrock test
			if(mpi_world_size() == 1)
			{
				Rosenbrock<Matrix, Vector> rosenbrock;
				Vector expected_rosenbrock = values(2, 1);

				auto subproblem = std::make_shared<SteihaugToint<Matrix, Vector> >();


				Vector x0 = values(2, 2.0);

				QuasiTrustRegion<Matrix, Vector> tr_solver(subproblem);
				auto cg = std::make_shared<ConjugateGradient<Matrix, Vector> >();
				tr_solver.set_linear_solver(cg);	

				auto hes_approx    = std::make_shared<SR1<Matrix, Vector> >();
				tr_solver.set_hessian_approximation_strategy(hes_approx);

				tr_solver.max_it(100); 
				tr_solver.verbose(true);
				tr_solver.solve(rosenbrock, x0);

				utopia_test_assert(approxeq(expected_rosenbrock, x0));
			}
		}


		SolverTest()
		: _n(10) { }

	private:
		int _n;

	};

	void runGenericSolversTest()
	{
		UTOPIA_UNIT_TEST_BEGIN("SolversTest");
#ifdef WITH_PETSC
		SolverTest<DMatrixd, DVectord, PetscScalar>().run();
#endif

// #ifdef WITH_BLAS
// 		SolverTest<Matrixd, Vectord, double>().run();
// #endif //WITH_BLAS

		UTOPIA_UNIT_TEST_END("SolversTest");
	}
}
