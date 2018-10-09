#include "utopia.hpp"
#include "utopia_SolverTest.hpp"
#include "utopia_assemble_laplacian_1D.hpp"
#include "test_problems/utopia_TestProblems.hpp"

namespace utopia
{
	template<class Matrix, class Vector>
	class QuasiNewtonTest 
	{
		public:
		
			void run_dense()
			{
				// UTOPIA_RUN_TEST(quasi_newton_test);
				// UTOPIA_RUN_TEST(Quasi_TR_test); 
			}

			void run_sparse()
			{
				UTOPIA_RUN_TEST(lbfgs_quasi_newton_test); 
			}			

			void quasi_newton_test()
			{
				// because dense matrices can not be sum-up in parallel
				if(mpi_world_size() > 1) return;
				
				Parameters params;
				params.atol(1e-9);
				params.rtol(1e-15);
				params.stol(1e-15);
				params.verbose(_verbose);
				
				auto lsolver = std::make_shared< ConjugateGradient<Matrix, Vector> >();
				auto hess_approx_BFGS   = std::make_shared<BFGS<Matrix, Vector> >();


				QuasiNewton<Matrix, Vector> nlsolver(hess_approx_BFGS, lsolver);
				nlsolver.set_parameters(params);

				auto line_search  = std::make_shared<utopia::Backtracking<Matrix, Vector> >();
				nlsolver.set_line_search_strategy(line_search);
				
				
				SimpleQuadraticFunction<Matrix, Vector> fun;
				
				Vector x = values(_n, 2.);
				Vector expected_1 = zeros(x.size());
				

				nlsolver.solve(fun, x);
				utopia_test_assert(approxeq(expected_1, x));
				
				TestFunctionND_1<Matrix, Vector> fun2(x.size().get(0));
				x = values(_n, 2.0);
				Vector expected_2 = values(x.size().get(0), 0.468919);

				nlsolver.solve(fun2, x);
				utopia_test_assert(approxeq(expected_2, x));

				Rosenbrock<Matrix, Vector> rosenbrock;
				Vector x0 = values(2, 0.5);
				nlsolver.solve(rosenbrock, x0);
				Vector expected_rosenbrock = values(2, 1.0);

				utopia_test_assert(approxeq(x0, expected_rosenbrock));
			}


			void Quasi_TR_test()
			{
				// rosenbrock test
				if(mpi_world_size() == 1)
				{
					Rosenbrock<Matrix, Vector> rosenbrock;
					Vector expected_rosenbrock = values(2, 1);

					auto subproblem = std::make_shared<SteihaugToint<Matrix, Vector> >();
					subproblem->set_preconditioner(std::make_shared<IdentityPreconditioner<Matrix, Vector> >());
					subproblem->atol(1e-10);

					Vector x0 = values(2, 2.0);

					QuasiTrustRegion<Matrix, Vector> tr_solver(subproblem);
					tr_solver.atol(1e-6); 
					tr_solver.rtol(1e-9);

					auto hes_approx   = std::make_shared<BFGS<Matrix, Vector> >();
					hes_approx->set_update_hessian(true); 

					tr_solver.set_hessian_approximation_strategy(hes_approx);

					tr_solver.max_it(100); 
					tr_solver.verbose(_verbose);
					tr_solver.delta0(1); 
					tr_solver.solve(rosenbrock, x0);

					utopia_test_assert(approxeq(expected_rosenbrock, x0));
				}
			}

			void lbfgs_quasi_newton_test()
			{
				std::cout<<"lbfgs_quasi_newton_test  \n"; 

				auto memory_size = 7; 

				Bratu1D<Matrix, Vector> fun(_n);
	    		Vector x = values(_n, -1.0);
	    		fun.apply_bc_to_initial_guess(x);

	    		auto linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector> >();
				auto hess_approx_BFGS   = std::make_shared<LBFGSB<Matrix,  Vector> >(memory_size, linear_solver);

				QuasiNewtonBound<Matrix, Vector> solver(hess_approx_BFGS, linear_solver);

				Vector lb   = local_values(local_size(x).get(0), -0.5);
				Vector ub   = local_values(local_size(x).get(0), 0.5);

				auto box = make_box_constaints(make_ref(lb), make_ref(ub));
	    		solver.set_box_constraints(box);				


	    		Vector grad = 0*x; 
	    		fun.gradient(x, grad); 
	    		Vector t, d; 

	    		hess_approx_BFGS->compute_breakpoints(grad, x, lb, ub, t); 


	    		// {
	    		// 	Write<Vector> w(t); 

	    		// 	auto r = range(t);
	    		// 	for(auto i=r.begin(); i < r.end(); i++)
	    		// 	{
	    		// 		if(i==0)
	    		// 			t.set(i,0); 
	    		// 		if(i==1)
	    		// 			t.set(i,-12); 	    				
	    		// 		if(i==2)
	    		// 			t.set(i,-5); 	    
	    		// 	}					    				
	    		// }


	    		auto t_current = 0.0; 
	    		hess_approx_BFGS->get_d_corresponding_to_ti(t, grad, d, t_current); 

	    		Vector sorted; 
	    		vec_unique_sort_serial(t, sorted, 3);  
	    		// disp(sorted); 



				
				// SimpleQuadraticFunction<DSMatrixd, DVectord> fun;

				// Parameters params;
				// params.atol(1e-9);
				// params.rtol(1e-15);
				// params.stol(1e-15);
				// params.verbose(_verbose);

				// const auto m = 3; 
				
				// auto linear_solver = std::make_shared<Factorization<DSMatrixd, DVectord>>();

				// auto hess_approx_BFGS   = std::make_shared<LBFGSB<DSMatrixd,  DVectord> >(memory_size, linear_solver);


		  // 		auto k = 15;

		  //       DVectord v = values(k, 999); 
		  //       DVectord y = values(k, 55); 
		  //       DVectord s = values(k, 1); 

		  //       hess_approx_BFGS->initialize(fun, v); 
		  //       hess_approx_BFGS->update(s, y); 

		  //       std::cout<<"---- solver end ---- \n"; 
					
			}

		QuasiNewtonTest()
		: _n(10), _verbose(true) { }
		
	private:
		int _n;
		bool _verbose; 
	};

	

	void runQuasiNewtonTest()
	{
		UTOPIA_UNIT_TEST_BEGIN("runQuasiNewtonTest");
		#ifdef WITH_PETSC
				QuasiNewtonTest<DMatrixd, DVectord>().run_dense();
				QuasiNewtonTest<DSMatrixd, DVectord>().run_sparse();
		#endif

		// #ifdef WITH_BLAS
		// 		QuasiNewtonTest<Matrixd, Vectord>().run_dense();
		// #endif //WITH_BLAS

		// #ifdef WITH_TRILINOS
		// 		QuasiNewtonTest<TSMatrixd, TVectord>().run_sparse();
		// #endif //WITH_TRILINOS				

		UTOPIA_UNIT_TEST_END("runQuasiNewtonTest");					
	}
}
