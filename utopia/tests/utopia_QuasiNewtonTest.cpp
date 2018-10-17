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
				// UTOPIA_RUN_TEST(lbfgs_quasi_newton_test); 
				UTOPIA_RUN_TEST(Quasi_TR_test_LBFGS); 
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


			void Quasi_TR_test_LBFGS()
			{
				auto memory_size = 7; 

				Bratu1D<Matrix, Vector> fun(_n);
	    		Vector x = values(_n, 0.0);
	    		fun.apply_bc_to_initial_guess(x);

	    		auto linear_solver = std::make_shared<GMRES<Matrix, Vector> >();
				auto hess_approx_BFGS   = std::make_shared<LBFGSB<Matrix, Matrix,  Vector> >(memory_size, linear_solver);


				auto subproblem = std::make_shared<SteihaugToint<Matrix, Vector> >();
				subproblem->set_preconditioner(std::make_shared<IdentityPreconditioner<Matrix, Vector> >());
				subproblem->atol(1e-10);

				QuasiTrustRegion<Matrix, Vector> tr_solver(subproblem);
				tr_solver.atol(1e-4); 
				tr_solver.rtol(1e-9);

				tr_solver.set_hessian_approximation_strategy(hess_approx_BFGS);

				tr_solver.max_it(100); 
				tr_solver.verbose(_verbose);
				tr_solver.delta0(1); 
				tr_solver.solve(fun, x);


				disp(x); 

			}			

			void lbfgs_quasi_newton_test()
			{
				auto memory_size = 5; 

				Bratu1D<Matrix, Vector> fun(_n);
	    		Vector x = values(_n, 0.0);
				Vector lb   = local_values(local_size(x).get(0), 5);
				Vector ub   = local_values(local_size(x).get(0), 10);			    		
	    		fun.apply_bc_to_initial_guess(x);


	    		auto linear_solver = std::make_shared<GMRES<Matrix, Vector> >();
	    		auto hess_approx_BFGS   = std::make_shared<LBFGSB<Matrix, Matrix,  Vector> >(memory_size, linear_solver);


	    		hess_approx_BFGS->initialize(fun, x); 


	    		Vector s = local_values(local_size(x).get(0), 2);
	    		Vector y = local_values(local_size(x).get(0), 5);

	    		hess_approx_BFGS->update(s, y); 
	    		hess_approx_BFGS->update(s, y); 
	    		hess_approx_BFGS->update(s, y); 
	    		// hess_approx_BFGS->update(s, y); 
	    		// hess_approx_BFGS->update(s, y); 
	    		// hess_approx_BFGS->update(s, y); 
	    		// hess_approx_BFGS->update(s, y); 	
	    		// hess_approx_BFGS->update(s, y); 
	    		// hess_approx_BFGS->update(s, y); 
	    		// hess_approx_BFGS->update(s, y); 	

	    		// auto L = hess_approx_BFGS->L_dots_; 
	    		// auto D = hess_approx_BFGS->d_elements_; 

	    		// std::cout<<"L_dots_: "<< L.size() << "  L_dots_[0].size(): "<< (hess_approx_BFGS->L_dots_[0]).size() << "  \n"; 


	    		// if(mpi_world_rank()==0)
	    		// {

		    	// 	for(auto i =0; i < memory_size; i++)
		    	// 	{
		    	// 		for(auto j =0; j < memory_size; j++)
		    	// 		{
		    	// 			std::cout<<"  "<< L[i][j]; 
		    	// 		}
		    	// 		std::cout<<"  \n"; 
		    	// 	}
		    	// }


		    	// std::cout<<"------------------- d ------------------- \n"; 

		    	// if(mpi_world_rank()==0)
	    		// {

		    	// 	for(auto i =0; i < memory_size; i++)
		    	// 	{

		    	// 		std::cout<<"  "<< D[i]; 
		    	// 	}
		    	// }

		    	// std::cout<<"   \n"; 



	    		// disp(hess_approx_BFGS->S_); 
	    		// disp(hess_approx_BFGS->Y_); 

	    		// disp(hess_approx_BFGS->W_); 





	   //  		{
	   //  			Write<Vector>  w(x); 

	   //  			auto r = range(x); 

	   //  			for(auto i=r.begin(); i != r.end(); ++i)
	   //  			{
	   //  				x.set(i, i); 
	   //  			}

	   //  		}

	   //  		disp(x);


	   //  		Vector feasible_set = local_zeros(local_size(x)); 

	   //  		{
	   //  			Write<Vector>  w(feasible_set); 

	   //  			auto r = range(feasible_set); 

	   //  			for(auto i=r.begin(); i != r.end(); ++i)
	   //  			{
	   //  				if(i==2 || i==3 || i==7 || i==9 || i==11 || i==14)
	   //  					feasible_set.set(i, 1.0); 
	   //  			}

	   //  		}


	   //  		// disp(feasible_set); 

				// ReducedPrimalMethod<Matrix, Vector> primal_method; 
				// auto local_size_feasible_set = primal_method.get_local_size_feasible_set(feasible_set); 
				// // std::cout<<"local_size: "<< local_size_feasible_set << "  \n"; 





				// Matrix M = values(size(feasible_set).get(0), 5, 1.0); 
				// disp(M); 

				// if(mpi_world_rank()==1)
				// 	std::cout<<"M: "<< local_size(M).get(0) <<"  M: "<< local_size(M).get(1) << "  \n"; 

				// auto global_feasible_set = sum(feasible_set); 

				// Matrix M; 

				// MPI_Comm comm = PetscObjectComm((PetscObject)raw_type(feasible_set));

    //     		MatCreate(comm, &raw_type(M));
    //     		MatSetFromOptions(raw_type(M));
    //     		MatSetType(raw_type(M), MATDENSE);
    //     		MatSetSizes(raw_type(M), local_size(feasible_set).get(0), 3, size(feasible_set).get(0), 3);
    //     		MatSetUp(raw_type(M));

    //     		MatSetOption(raw_type(M), MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE);
    //     		MatSetOption(raw_type(M), MAT_NO_OFF_PROC_ENTRIES,     PETSC_FALSE);


    //     		{
    //     			Write<Matrix>  mr(M); 

    //     			auto r_range = row_range(M); 
    //     			auto c_range = col_range(M); 


    //     			for(auto i=r_range.begin(); i != r_range.end(); ++i)
    //     			{
    //     				for(auto j=c_range.begin(); j != c_range.end(); ++j)
	   //      			{
	   //      				M.set(i,j, j); 
	   //      			}
    //     			}
    //     		}

    //     		disp(M); 

				// if(mpi_world_rank()==0)
				// 	std::cout<<"M: "<< local_size(M).get(0) <<"  M: "<< local_size(M).get(1) << "  \n";         		

				// Matrix M_T = transpose(M); 
				// // disp(M_T); 

				// if(mpi_world_rank()==0)
				// 	std::cout<<"MT: "<< local_size(M_T).get(0)<<"  MT: "<< local_size(M_T).get(1) << "  \n"; 


				// Matrix M_reduced  = local_values(local_size(M_T).get(0), local_size_feasible_set, 0.0); 

				// disp(M_reduced); 

				// primal_method.build_reduced_matrix(M_T, feasible_set, M_reduced); 

				// disp(M_reduced); 



	   			// auto linear_solver = std::make_shared<GMRES<Matrix, Vector> >();
				// auto hess_approx_BFGS   = std::make_shared<LBFGSB<Matrix, Matrix,  Vector> >(memory_size, linear_solver);


				// auto subproblem = std::make_shared<SteihaugToint<Matrix, Vector> >();
				// subproblem->set_preconditioner(std::make_shared<IdentityPreconditioner<Matrix, Vector> >());
				// subproblem->atol(1e-10);

				// QuasiTrustRegion<Matrix, Vector> tr_solver(subproblem);
				// tr_solver.atol(1e-4); 
				// tr_solver.rtol(1e-9);

				// tr_solver.set_hessian_approximation_strategy(hess_approx_BFGS);

				// tr_solver.max_it(100); 
				// tr_solver.verbose(_verbose);
				// tr_solver.delta0(1); 
				// tr_solver.solve(fun, x);
				// disp(x); 








				// QuasiNewtonBound<Matrix, Vector> solver(hess_approx_BFGS, linear_solver);

				// auto line_search  = std::make_shared<utopia::Backtracking<Matrix, Vector> >();
				// solver.set_line_search_strategy(line_search);
				// solver.max_it(10); 


				// Vector lb   = local_values(local_size(x).get(0), -0.01);
				// Vector ub   = local_values(local_size(x).get(0), 0.01);				

				// auto box = make_box_constaints(make_ref(lb), make_ref(ub));
	   //  		solver.set_box_constraints(box);				


	   //  		solver.verbose(true); 
	   //  		solver.solve(fun, x);
	   //  		disp(x); 	



	   //  		DMatrixd M; 
	   //  		M=zeros(10,10); 
	   //  		assemble_symmetric_laplacian_1D(M, true); 
    //             disp(M); 

				
				// auto gmres = std::make_shared<GMRES<DMatrixd, DVectord> >();
				// MatLinearSolver<DMatrixd, DMatrixd, DVectord> mat_solver(gmres); 

				// DMatrixd M_inv; 
				// mat_solver.get_inverse(M, M_inv); 

				// disp(M_inv); 


				// DMatrixd MM_iinnvv = inv(M); 

				// disp(MM_iinnvv); 



			   	//  		Vector result; 
				//  		hess_approx_BFGS->apply_H(x, result); 
				//  		disp(result); 
				// GeneralizedCauchyPoint<Matrix, Vector> cp; 
				// cp.print(); 


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
				// QuasiNewtonTest<DMatrixd, DVectord>().run_dense();
			
			QuasiNewtonTest<DSMatrixd, DVectord>().run_sparse();
			// QuasiNewtonTest<DMatrixd, DVectord>().run_sparse();
		#endif

		#ifdef WITH_BLAS
				// QuasiNewtonTest<Matrixd, Vectord>().run_dense();
		#endif //WITH_BLAS

		// #ifdef WITH_TRILINOS
		// 		QuasiNewtonTest<TSMatrixd, TVectord>().run_sparse();
		// #endif //WITH_TRILINOS				

		UTOPIA_UNIT_TEST_END("runQuasiNewtonTest");					
	}
}
