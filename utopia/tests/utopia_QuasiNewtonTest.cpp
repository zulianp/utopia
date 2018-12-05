#include "utopia.hpp"
#include "utopia_SolverTest.hpp"
#include "utopia_assemble_laplacian_1D.hpp"
#include "test_problems/utopia_TestProblems.hpp"
#include "test_problems/utopia_BratuMultilevelTestProblem.hpp"

namespace utopia
{
	template<class Matrix, class Vector, class ApproxType>
	class QuasiNewtonTest 
	{

		public:
			static void print_backend_info()
			{
				if(Utopia::instance().verbose() && mpi_world_rank() == 0) {
					std::cout << "\nBackend: " << backend_info(Vector()).get_name() << std::endl;
				}
			}
		
			void run_dense()
			{	
				UTOPIA_RUN_TEST(quasi_newton_test);
				UTOPIA_RUN_TEST(Quasi_TR_test); 
			}

			void run_sparse()
			{
				UTOPIA_RUN_TEST(Quasi_TR_test_sparse); 
				UTOPIA_RUN_TEST(quasi_newton_test_sparse); 
				UTOPIA_RUN_TEST(QuasiTR_constraint_GCP_test);
				UTOPIA_RUN_TEST(Quasi_TR_Gradient_projection_active_set_test); 
				UTOPIA_RUN_TEST(QuasiNewtonBoundTest); 
				
				UTOPIA_RUN_TEST(TR_constraint_GCP_test);
				UTOPIA_RUN_TEST(Gradient_projection_active_set_test);				
			}	

			void run_multilevel()
			{
				UTOPIA_RUN_TEST(Quasi_RMTR_test); 
				UTOPIA_RUN_TEST(Quasi_RMTR_inf_bound_test); 
			}	

			void quasi_newton_test()
			{
				// because dense matrices can not be sum-up in parallel
				if(mpi_world_size() > 1) return;
				
				Parameters params;
				params.atol(1e-5);
				params.rtol(1e-15);
				params.stol(1e-15);
				params.verbose(_verbose);
				
				auto hessian_approx   = std::make_shared<ApproxType >();
				auto lsolver = std::make_shared<EmptyPrecondMatrixFreeLinearSolver<Vector> >(); 

				auto precond = hessian_approx->build_Hinv_precond(); 
		        lsolver->set_preconditioner(precond); 

				QuasiNewton<Vector> nlsolver(hessian_approx, lsolver);
				nlsolver.set_parameters(params);

				auto line_search  = std::make_shared<utopia::Backtracking<Vector> >();
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

					auto subproblem = std::make_shared<SteihaugToint<Matrix, Vector, HOMEMADE> >();
					subproblem->set_preconditioner(std::make_shared<IdentityPreconditioner<Matrix, Vector> >());
					subproblem->atol(1e-10);

					Vector x0 = values(2, 2.0);

					auto hes_approx   = std::make_shared<ApproxType >();
					hes_approx->set_update_hessian(true); 					

					QuasiTrustRegion<Vector> tr_solver(hes_approx, subproblem);
					tr_solver.atol(1e-6); 
					tr_solver.rtol(1e-9);

					tr_solver.max_it(1000); 
					tr_solver.verbose(_verbose);
					tr_solver.delta0(1); 
					tr_solver.solve(rosenbrock, x0);

					utopia_test_assert(approxeq(expected_rosenbrock, x0));
				}
			}			

			void Quasi_TR_test_sparse()
			{
				auto memory_size = 7; 

				Bratu1D<Matrix, Vector> fun(_n);
	    		Vector x = values(_n, 1.0);
	    		fun.apply_bc_to_initial_guess(x);

				auto hess_approx = std::make_shared<ApproxType >(memory_size);
				auto subproblem = std::make_shared<SteihaugToint<Matrix, Vector, HOMEMADE> >();
				subproblem->set_preconditioner(std::make_shared<IdentityPreconditioner<Matrix, Vector> >());
				// subproblem->verbose(true);

				QuasiTrustRegion<Vector> tr_solver(hess_approx, subproblem);
				tr_solver.atol(1e-5); 
				tr_solver.rtol(1e-9);

				tr_solver.max_it(2000); 
				tr_solver.verbose(_verbose);
				tr_solver.delta0(1); 
				tr_solver.solve(fun, x);
			}					

			void quasi_newton_test_sparse()
			{	
				SizeType memory_size = 5; 			

				Parameters params;
				params.atol(1e-6);
				params.rtol(1e-15);
				params.stol(1e-15);
				params.max_it(1000); 
				params.verbose(_verbose);
				
	    		auto hess_approx   = std::make_shared<ApproxType >(memory_size);				
	    		auto lsolver = std::make_shared<EmptyPrecondMatrixFreeLinearSolver<Vector> >(); 
	    		
	    		auto precond = hess_approx->build_Hinv_precond(); 
		        lsolver->set_preconditioner(precond); 


				QuasiNewton<Vector> nlsolver(hess_approx, lsolver);
				nlsolver.set_parameters(params);

				auto line_search  = std::make_shared<utopia::Backtracking<Vector> >();
				nlsolver.set_line_search_strategy(line_search);
				
				SimpleQuadraticFunction<Matrix, Vector> fun;
				
				Vector x = values(_n, 2.);
				Vector expected_1 = zeros(x.size());
				
				nlsolver.solve(fun, x);
				utopia_test_assert(approxeq(expected_1, x));	


				Bratu1D<Matrix, Vector> fun2(_n);
	    		Vector x2 = values(_n, 1.0); 		
	    		fun2.apply_bc_to_initial_guess(x2);
	    		nlsolver.solve(fun2, x2);

			}
		

		    void TR_constraint_GCP_test()
		    {
		    	Bratu1D<Matrix, Vector> fun(_n);
		    	Vector x = values(_n, 1.0);
		    	fun.apply_bc_to_initial_guess(x);

		    	DVectord ub, lb;
		    	fun.generate_constraints(lb, ub);
		    	auto box = make_box_constaints(make_ref(lb), make_ref(ub));

		    	Parameters params;
				params.atol(1e-6);
				params.rtol(1e-10);
				params.stol(1e-10);
				params.max_it(1000); 
				params.verbose(_verbose);

		        auto qp_solver = std::make_shared<GeneralizedCauchyPoint<Matrix, Vector> >();

		        TrustRegionVariableBound<Matrix, Vector>  tr_solver(qp_solver);
		        tr_solver.set_box_constraints(box);
				tr_solver.set_parameters(params);
				tr_solver.solve(fun, x);
		    }

		    void QuasiTR_constraint_GCP_test()
		    {
		    	Bratu1D<Matrix, Vector> fun(_n);
		    	Vector x = values(_n, 1.0);
		    	fun.apply_bc_to_initial_guess(x);

				Vector lb   = local_values(local_size(x).get(0), -0.01);
				Vector ub   = local_values(local_size(x).get(0), 0.01);		
		    	auto box = make_box_constaints(make_ref(lb), make_ref(ub));

		    	SizeType memory_size = 5; 
		    	Parameters params;
				params.atol(1e-6);
				params.rtol(1e-10);
				params.stol(1e-10);
				params.verbose(_verbose);
				params.max_it(1000);
				params.delta0(1);

				auto hess_approx   = std::make_shared<ApproxType >(memory_size);	
		        auto qp_solver = std::make_shared<GeneralizedCauchyPoint<Matrix, Vector> >();
		        qp_solver->set_memory_size(10); 

		        QuasiTrustRegionVariableBound<Vector>  tr_solver(hess_approx, qp_solver);
		        tr_solver.set_box_constraints(box);
				tr_solver.set_parameters(params);
				tr_solver.solve(fun, x);

		    }


		    void Gradient_projection_active_set_test()
		    {
		    	Bratu1D<Matrix, Vector> fun(_n);
		    	Vector x = values(_n, 0.0);
		    	fun.apply_bc_to_initial_guess(x);

				Vector lb   = local_values(local_size(x).get(0), -0.01);
				Vector ub   = local_values(local_size(x).get(0), 0.01);		
		    	auto box = make_box_constaints(make_ref(lb), make_ref(ub));

		    	Parameters params;
				params.atol(1e-6);
				params.rtol(1e-10);
				params.stol(1e-10);
				params.verbose(_verbose);
				params.max_it(1000);
				params.delta0(1); 

		        auto qp_solver = std::make_shared<ProjectedGradientActiveSet<Matrix, Vector> >();
		        qp_solver->verbose(false); 
		        qp_solver->atol(1e-12); 

		        auto lsolver = std::make_shared<LUDecomposition<Matrix, Vector> >();
		        qp_solver->set_linear_solver(lsolver); 


		        TrustRegionVariableBound<Matrix, Vector>  tr_solver(qp_solver);
		        tr_solver.set_box_constraints(box);
				tr_solver.set_parameters(params);
				tr_solver.solve(fun, x);

		    }



		    void Quasi_TR_Gradient_projection_active_set_test()
		    {
		    	SizeType memory_size = 5; 

		    	Bratu1D<Matrix, Vector> fun(_n);
		    	Vector x = values(_n, 0.0);
		    	fun.apply_bc_to_initial_guess(x);

				Vector lb   = local_values(local_size(x).get(0), -0.01);
				Vector ub   = local_values(local_size(x).get(0), 0.01);		
		    	auto box = make_box_constaints(make_ref(lb), make_ref(ub));

		    	Parameters params;
				params.atol(1e-6);
				params.rtol(1e-10);
				params.stol(1e-10);
				params.verbose(_verbose);
				params.max_it(20); 
				params.delta0(1); 

				auto hess_approx   = std::make_shared<ApproxType >(memory_size);	
		        auto qp_solver = std::make_shared<ProjectedGradientActiveSet<Matrix, Vector> >();

				auto precond = hess_approx->build_Hinv_precond(); 
		        qp_solver->set_preconditioner(precond); 

		        QuasiTrustRegionVariableBound<Vector>  tr_solver(hess_approx, qp_solver);
		        tr_solver.set_box_constraints(box);
				tr_solver.set_parameters(params);
				tr_solver.solve(fun, x);
		    }

			void QuasiNewtonBoundTest()
			{
				auto memory_size = 5; 

				Bratu1D<Matrix, Vector> fun(_n);
	    		Vector x = values(_n, 0.0);
				Vector lb   = local_values(local_size(x).get(0), -0.01);
				Vector ub   = local_values(local_size(x).get(0), 0.01);		    		
	    		fun.apply_bc_to_initial_guess(x);

				auto hess_approx   = std::make_shared<ApproxType >(memory_size);	
		        auto qp_solver = std::make_shared<ProjectedGradientActiveSet<Matrix, Vector> >();
		        
				QuasiNewtonBound<Vector> solver(hess_approx, qp_solver);

				auto line_search  = std::make_shared<utopia::Backtracking<Vector> >();
				solver.set_line_search_strategy(line_search);
				solver.max_it(50); 	
				solver.stol(1e-12); 	
				solver.atol(1e-6); 		

				auto box = make_box_constaints(make_ref(lb), make_ref(ub));
	    		solver.set_box_constraints(box);				

	    		solver.verbose(_verbose); 
	    		solver.solve(fun, x);
			}

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		    void Quasi_RMTR_test()
		    {	
		    	const SizeType n_levels = 3; 
		    	BratuMultilevelTestProblem<Matrix, Vector> problem(n_levels, true, _verbose); 

		    	auto rmtr = std::make_shared<QuasiRMTR<Matrix, Vector, FIRST_ORDER>  >(n_levels);
		    
		    	// intial guess
		        Vector x = values(problem.n_dofs[problem.n_levels -1 ], 0.0);
		    	std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector> > >  level_functions(problem.n_levels);

		    	for(auto l=0; l < problem.n_levels; l++)
		    	{
			    	auto fun = std::make_shared<Bratu1D<Matrix, Vector> >(problem.n_dofs[l]);
			    	level_functions[l] = fun;

			    	// making sure that fine level IG is feasible
			    	if(l+1 == problem.n_levels)
			    		fun->apply_bc_to_initial_guess(x);
			    }

		        const SizeType memory_size = 5; 
		        std::vector<std::shared_ptr<HessianApproximation<Vector> > > hess_approxs(problem.n_levels);
		    	for(auto l=0; l < problem.n_levels; l++)
		    	{
		    		auto hes_approx   = std::make_shared<ApproxType >(memory_size);
		    		hess_approxs[l] = hes_approx; 
			    }

				rmtr->set_hessian_approximation_strategies(hess_approxs);

			    std::vector<std::shared_ptr<utopia::MatrixFreeTRSubproblem<Vector> > > subproblems(problem.n_levels);
		    	for(auto l=0; l < problem.n_levels; l++)
		    	{
		    		auto tr_strategy = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();

		    		auto precond = hess_approxs[l]->build_Hinv_precond(); 
		        	tr_strategy->set_preconditioner(precond); 

		    		subproblems[l] = tr_strategy; 
			    }			 

		        rmtr->set_tr_strategies(subproblems); 			     
		        rmtr->set_transfer_operators(problem.prolongations, problem.restrictions);

		        rmtr->max_it(50);
		        rmtr->max_coarse_it(3);
		        rmtr->max_smoothing_it(3);
		        rmtr->delta0(1);
		        rmtr->atol(1e-4);
		        rmtr->rtol(1e-10);
		        rmtr->set_grad_smoothess_termination(0.000001);
		        rmtr->set_eps_grad_termination(1e-7);

				rmtr->verbose(problem.verbose);
				// rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
				rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);
		        rmtr->set_functions(level_functions);


		        rmtr->solve(x);
		    }



		void Quasi_RMTR_inf_bound_test()
	    {
	    	const SizeType n_levels = 3; 
		    BratuMultilevelTestProblem<Matrix, Vector> problem(n_levels, true, _verbose); 

        	auto rmtr = std::make_shared<QuasiRMTR_inf<Matrix, Vector, FIRST_ORDER>  >(n_levels);

	    	// intial guess
	        Vector x = values(problem.n_dofs[problem.n_levels -1 ], 0.0);

	        // upper, lower bound...
	        Vector ub, lb;
	        std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector> > >  level_functions(problem.n_levels);
	    	for(auto l=0; l < problem.n_levels; l++)
	    	{
		    	auto fun = std::make_shared<Bratu1D<Matrix, Vector> >(problem.n_dofs[l]);
		    	level_functions[l] = fun;

		    	// making sure that fine level IG is feasible
		    	if(l+1 == problem.n_levels)
		    	{
		    		fun->apply_bc_to_initial_guess(x);
		    		fun->generate_constraints(lb, ub, -10, 0.01);
		    	}
		    }

	        const SizeType memory_size = 5; 
	        std::vector<std::shared_ptr<HessianApproximation<Vector> > > hess_approxs(problem.n_levels);
	    	for(auto l=0; l < problem.n_levels; l++)
	    	{
	    		auto hes_approx   = std::make_shared<ApproxType >(memory_size);
	    		hess_approxs[l] = hes_approx; 
		    }

			rmtr->set_hessian_approximation_strategies(hess_approxs);


		    std::vector<std::shared_ptr<utopia::MatrixFreeQPSolver<Vector> > > subproblems(problem.n_levels);
	    	for(auto l=0; l < problem.n_levels; l++)
	    	{
	    		auto tr_strategy = std::make_shared<utopia::ProjectedGradientActiveSet<Matrix, Vector> >();
	    		auto precond = hess_approxs[l]->build_Hinv_precond(); 
	        	tr_strategy->set_preconditioner(precond); 

	    		subproblems[l] = tr_strategy; 
		    }	

	        rmtr->set_tr_strategies(subproblems); 
	        rmtr->set_transfer_operators(problem.prolongations, problem.restrictions);

	        rmtr->max_it(30);
	        rmtr->max_coarse_it(3);
	        rmtr->max_smoothing_it(3);
	        rmtr->delta0(1);
	        rmtr->atol(1e-4);
	        rmtr->rtol(1e-10);
	        rmtr->set_grad_smoothess_termination(0.000001);
	        rmtr->set_eps_grad_termination(1e-7);

			rmtr->verbose(problem.verbose);
			// rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
			rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);

	        rmtr->set_functions(level_functions);

	       	auto box = make_box_constaints(make_ref(lb), make_ref(ub));
	    	rmtr->set_box_constraints(box);

	        rmtr->solve(x);

	    }




		QuasiNewtonTest()
		: _n(10), _verbose(false) { }
		
	private:
		int _n;
		bool _verbose; 
	};

	

	void runQuasiNewtonTest()
	{
		UTOPIA_UNIT_TEST_BEGIN("runQuasiNewtonTest");
		#ifdef WITH_PETSC
			QuasiNewtonTest<DMatrixd, DVectord, BFGS<DMatrixd, DVectord> >().print_backend_info();
			QuasiNewtonTest<DMatrixd, DVectord, BFGS<DMatrixd, DVectord> >().run_dense();
			
			QuasiNewtonTest<DSMatrixd, DVectord, LBFGS<DVectord> >().run_sparse();
			QuasiNewtonTest<DSMatrixd, DVectord, LSR1<DVectord> >().run_sparse();

			QuasiNewtonTest<DSMatrixd, DVectord, LBFGS<DVectord> >().run_multilevel();
		#endif

		#ifdef WITH_BLAS
				QuasiNewtonTest<Matrixd, Vectord, BFGS<Matrixd, Vectord> >().print_backend_info();
				QuasiNewtonTest<Matrixd, Vectord, BFGS<Matrixd, Vectord> >().run_dense();
		#endif //WITH_BLAS

		// #ifdef WITH_TRILINOS
				// QuasiNewtonTest<TSMatrixd, TVectord>().print_backend_info();
		// 		QuasiNewtonTest<TSMatrixd, TVectord>().run_sparse();
		// #endif //WITH_TRILINOS				

		UTOPIA_UNIT_TEST_END("runQuasiNewtonTest");					
	}
}
