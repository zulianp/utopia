/*
* @Author: kopanicakova
* @Date:   2018-02-06 17:47:26
* @Last Modified by:   kopanicakova
* @Last Modified time: 2018-04-10 11:02:22
*/
#include "utopia.hpp"
#include "utopia_SolverTest.hpp"
#include "test_problems/utopia_TestProblems.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"

namespace utopia
{

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
			UTOPIA_RUN_TEST(ngs_test);
			UTOPIA_RUN_TEST(newton_cg_test);
			UTOPIA_RUN_TEST(solver_from_params_test);
			UTOPIA_RUN_TEST(tr_test);
			UTOPIA_RUN_TEST(ls_test);
			UTOPIA_RUN_TEST(nl_solve_test);
			
		}
		
		class EmptyLSFun : public LeastSquaresFunction<Matrix, Vector> {
		public:
			bool residual(const Vector &/*point*/, Vector &/*result*/) const {}
			bool jacobian(const Vector &/*x*/, Matrix &/*hessian*/) const {}
			bool value(const Vector &, Scalar &val) const {}
			bool update(const Vector &) { }
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
			assert(approxeq(expected, actual));
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
			assert(approxeq(expected, actual));
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
				// assert(approxeq(expected, x));
				
				
				x = values(10, 2);
				params.trust_region_alg(STEIHAUG_TOINT_TAG);
				trust_region_solve(fun2, x, params);
				assert(approxeq(expected, x));
				
				
				// x = values(10, 2);
				// params.trust_region_alg(STEIHAUG_TOINT_TAG);
				
				// auto strategy_tr = std::make_shared<utopia::SteihaugToint<Matrix, Vector> >();
				// auto precond = std::make_shared< InvDiagPreconditioner<Matrix, Vector> >();
				// strategy_tr->set_preconditioner(precond);
				// strategy_tr->verbose(true);
				// TrustRegion<Matrix, Vector> tr_solver(strategy_tr);
				
				// tr_solver.verbose(true);
				// tr_solver.solve(fun2, x);
				// assert(approxeq(expected, x));
				
				
				
				x = values(10, 2);
				params.trust_region_alg(CAUCHYPOINT_TAG);
				trust_region_solve(fun2, x, params);
				assert(approxeq(expected, x));
				
				x = values(10, 2);
				trust_region_solve(fun2, x, params);
				assert(approxeq(expected, x));
				
				
				Vector expected_rosenbrock = values(2, 1);
				
				Rosenbrock<Matrix, Vector> rosenbrock;
				// params.trust_region_alg(DOGLEG_TAG);
				Vector x0 = values(2, 2.0);
				// trust_region_solve(rosenbrock, x0, params);
				// assert(approxeq(expected_rosenbrock, x0));
				
				x0 = values(2, 2.0);
				params.trust_region_alg(STEIHAUG_TOINT_TAG);
				trust_region_solve(rosenbrock, x0, params);
				assert(approxeq(expected_rosenbrock, x0));
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
				
				assert(approxeq(expected_woods, x_w1));
				assert(approxeq(expected_woods, x_w2));
				
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
				assert(approxeq(expected_rosenbrock, x01));
				assert(approxeq(expected_rosenbrock, x02));
				
				x01 = values(2, 2.0);
				params.verbose(false);
				params.linear_solver_verbose(false);
				line_search_solve(rosenbrock_fun, x01, params);
			}
		}

		void ngs_test()
		{			
			const SizeType n = 30;

			Matrix m = zeros(n, n);
			assemble_laplacian_1D(n, m);
			{
				Range r = row_range(m);
				Write<Matrix> w(m);
				if(r.begin() == 0) {
					m.set(0, 0, 1.);
					m.set(0, 1, 0);
				}

				if(r.end() == n) {
					m.set(n-1, n-1, 1.);
					m.set(n-1, n-2, 0);
				}
			}

			Vector rhs = values(n, 1.);
			{ 
			    //Creating test vector (alternative way see [assemble vector alternative], which might be easier for beginners)
				Range r = range(rhs);
				Write<Vector> w(rhs);

				if(r.begin() == 0) {
					rhs.set(0, 0);
				}

				if(r.end() == n) {
					rhs.set(n-1, 0.);
				}
			}

			Vector upper_bound = values(n, 100.0);
			Vector solution    = zeros(n);

			ProjectedGaussSeidel<Matrix, Vector> pgs;
			pgs.max_it(n*2);
			// pgs.verbose(true);
			pgs.set_box_constraints(make_upper_bound_constraints(make_ref(upper_bound)));

			Chrono c;
			c.start();
			pgs.solve(m, rhs, solution);
			c.stop();

			// std::cout << c << std::endl;
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
		
#ifdef WITH_BLAS
		SolverTest<Matrixd, Vectord, double>().run();
#endif //WITH_BLAS
		
		UTOPIA_UNIT_TEST_END("SolversTest");
	}


    void runSolversTest()
    {
    	runGenericSolversTest(); 
    	runPetscNonlinearSolversTest(); 
    	runPetscLinearSolversTest(); 
    	runPetscSlepcSolversTest(); 
    }


}
