#include "utopia_Testing.hpp"
#include "utopia_Chrono.hpp"
#include "utopia_MPI.hpp"
#include "utopia.hpp"
#include "utopia_Base.hpp"
#include "utopia_Benchmark.hpp"
#include "utopia_LargeScaleIncludes.hpp"
#include "utopia_UnconstrainedBenchmark.hpp"
#include "utopia_ConstrainedBenchmark.hpp"
#include <string>
#include <cassert>

namespace utopia 
{

	template<class Matrix, class Vector>
	class QuasiTRTest : public Benchmark 
	{
	public:
		DEF_UTOPIA_SCALAR(Vector);

		virtual std::string name() override
		{
			return "QuasiTRTest benchmark.";
		}

		QuasiTRTest(const SizeType & n = 10, const bool verbose = false): n_(n), verbose_(verbose)
		{
			if(mpi_world_size()==1)
			{
				test_functions_parallel_.resize(4);
				test_functions_parallel_[0] = std::make_shared<Poisson1D<Matrix, Vector> >(n_);
				test_functions_parallel_[1] = std::make_shared<Bratu1D<Matrix, Vector> >(n_);
				test_functions_parallel_[2] = std::make_shared<Morebv1D<Matrix, Vector> >(n_*mpi_world_size());
				test_functions_parallel_[3] = std::make_shared<NonEllipse2D<Matrix, Vector> >(n_*mpi_world_size());				

				test_functions_.resize(4);
				test_functions_[0] = std::make_shared<Rosenbrock01<Matrix, Vector> >();
				test_functions_[1] = std::make_shared<QPTestFunction_2D<Matrix, Vector> >();
				test_functions_[2] = std::make_shared<Biggs18<Matrix, Vector> >();
				test_functions_[3] = std::make_shared<VariablyDim25<Matrix, Vector> >();

				test_functions_constrained_.resize(2); 
				test_functions_constrained_[0] = std::make_shared<Beale05Constrained<Matrix, Vector> >();
				// test_functions_constrained_[1] = std::make_shared<Watson20Constrained<Matrix, Vector> >();
				// test_functions_constrained_[1] = std::make_shared<VariablyDim25Constrained<Matrix, Vector> >();
				// test_functions_constrained_[3] = std::make_shared<PenaltyII24Constrained<Matrix, Vector> >();
				test_functions_constrained_[1] = std::make_shared<Rosenbrock21Constrained<Matrix, Vector> >();
			}
			else
			{
				test_functions_parallel_.resize(4);
				test_functions_parallel_[0] = std::make_shared<Poisson1D<Matrix, Vector> >(n_*mpi_world_size());
				test_functions_parallel_[1] = std::make_shared<Bratu1D<Matrix, Vector> >(n_*mpi_world_size());
				test_functions_parallel_[2] = std::make_shared<Morebv1D<Matrix, Vector> >(n_*mpi_world_size());

				// Only petsc
				test_functions_parallel_[3] = std::make_shared<NonEllipse2D<Matrix, Vector> >(n_*mpi_world_size());
			}
		}

		~QuasiTRTest()
		{
			test_functions_.clear(); 
		}


		void initialize() override
		{


			// this->register_experiment("QuasiNewton_InvApprox_LBFGS",
			// 	[this]() {
	  //               auto hessian_approx   = std::make_shared<LBFGS<Vector> >(7);
	  //               auto lsolver = std::make_shared<EmptyPrecondMatrixFreeLinearSolver<Vector> >();

	  //               auto precond = hessian_approx->build_Hinv_precond();
	  //               lsolver->set_preconditioner(precond);

	  //               QuasiNewton<Vector> solver(hessian_approx, lsolver);

	  //               auto line_search  = std::make_shared<utopia::Backtracking<Vector> >();
	  //               solver.set_line_search_strategy(line_search);

		 //            run_tr(this->test_functions_parallel_, solver, "QuasiNewton_InvApprox_LBFGS", this->verbose_);
		 //            run_tr(this->test_functions_, solver, "QuasiNewton_InvApprox_LBFGS", this->verbose_);
			// 	}
			// );

			// this->register_experiment("QuasiNewton_InvApprox_BFGS",
			// 	[this]() {
	  //               auto hessian_approx   = std::make_shared<BFGS<Matrix, Vector> >();
	  //               auto lsolver = std::make_shared<EmptyPrecondMatrixFreeLinearSolver<Vector> >();

	  //               auto precond = hessian_approx->build_Hinv_precond();
	  //               lsolver->set_preconditioner(precond);

	  //               QuasiNewton<Vector> solver(hessian_approx, lsolver);

	  //               auto line_search  = std::make_shared<utopia::Backtracking<Vector> >();
	  //               solver.set_line_search_strategy(line_search);

		 //            run_tr(this->test_functions_parallel_, solver, "QuasiNewton_InvApprox_BFGS", this->verbose_);
		 //            run_tr(this->test_functions_, solver, "QuasiNewton_InvApprox_BFGS", this->verbose_);
			// 	}
			// );			


			// this->register_experiment("QuasiTR_LBFGS",
			// 	[this]() {
	  //               auto hessian_approx   = std::make_shared<LBFGS<Vector> >(7);
   //                  auto subproblem = std::make_shared<SteihaugToint<Matrix, Vector, HOMEMADE> >();

   //                  subproblem->set_preconditioner(std::make_shared<IdentityPreconditioner<Vector> >());
   //                  // subproblem->atol(1e-10);

   //                  QuasiTrustRegion<Vector> solver(hessian_approx, subproblem);

		 //            run_tr(this->test_functions_parallel_, solver, "QuasiTR_LBFGS", this->verbose_);
		 //            run_tr(this->test_functions_, solver, "QuasiTR_LBFGS", this->verbose_);
			// 	}
			// );


			// this->register_experiment("QuasiTR_LSR1",
			// 	[this]() {
	  //               auto hessian_approx   = std::make_shared<LSR1<Vector> >(7);
   //                  auto subproblem = std::make_shared<SteihaugToint<Matrix, Vector, HOMEMADE> >();

   //                  subproblem->set_preconditioner(std::make_shared<IdentityPreconditioner<Vector> >());
   //                  // subproblem->atol(1e-10);

   //                  QuasiTrustRegion<Vector> solver(hessian_approx, subproblem);

		 //            run_tr(this->test_functions_parallel_, solver, "QuasiTR_LSR1", this->verbose_);
		 //            run_tr(this->test_functions_, solver, "QuasiTR_LSR1", this->verbose_);
			// 	}
			// );
					


			// this->register_experiment("QuasiNewton_Bound_LBFGS",
			// 	[this]() {
	  //               auto hessian_approx   = std::make_shared<LBFGS<Vector> >(7);
	  //               auto qp_solver = std::make_shared<MPGRP<Matrix, Vector> >();

	  //               QuasiNewtonBound<Vector> solver(hessian_approx, qp_solver);

	  //               auto line_search  = std::make_shared<utopia::Backtracking<Vector> >();
	  //               solver.set_line_search_strategy(line_search);

		 //            run_tr(this->test_functions_constrained_, solver, "QuasiNewton_Bound_LBFGS", this->verbose_);
			// 		run_tr(this->test_functions_, solver, "QuasiNewton_Bound_LBFGS", this->verbose_);
			// 		//run_tr(this->test_functions_parallel_, solver, "QuasiNewton_Bound_LBFGS", this->verbose_);
			// 	}
			// );


			this->register_experiment("QuasiTR_Bound_LBFGS",
				[this]() {
	                auto hessian_approx   = std::make_shared<LBFGS<Vector> >(7);
	                auto qp_solver = std::make_shared<MPGRP<Matrix, Vector> >();
	                qp_solver->max_it(1000);
	                qp_solver->atol(1e-12);
	                // qp_solver->max_it(10);
	                // qp_solver->verbose(true);

	                // QuasiNewtonBound<Vector> solver(hessian_approx, qp_solver);
	                QuasiTrustRegionVariableBound<Vector>  solver(hessian_approx, qp_solver);
	                // solver.delta0(1e-5);


		            // run_tr(this->test_functions_constrained_, solver, "QuasiTR_Bound_LBFGS", this->verbose_);
					// run_tr(this->test_functions_, solver, "QuasiTR_Bound_LBFGS", this->verbose_);
					run_tr(this->test_functions_parallel_, solver, "QuasiNewton_Bound_LBFGS", this->verbose_);
				}
			);
	
		}

	private:

		template<class Fun, class NonlinearSolver>
		static void run_tr(std::vector<std::shared_ptr<Fun> > & test_functions, NonlinearSolver &solver, const std::string & solv_name,  const bool & exp_verbose = false) 
		{
			InputParameters in;
			in.set("atol", 1e-5);
			in.set("rtol", 1e-11);
			in.set("stol", 1e-14);
			in.set("stol", 1e-14);
			in.set("delta_min", 1e-13); 
			in.set("max-it", 100); 
			in.set("verbose", true);
			solver.read(in); 

			if(exp_verbose && mpi_world_rank()==0)
			{
				std::cout<<"--------------------------------------------------------- \n";
				std::cout<<"				" << solv_name << "				\n";
				std::cout<<"--------------------------------------------------------- \n";
			}

	    	for(size_t i =0; i < test_functions.size(); i++)
	    	{
				Vector x_init = test_functions[i]->initial_guess(); 

				if(dynamic_cast<ConstrainedTestFunction<Matrix, Vector> *>(test_functions[i].get()))
				{
					ConstrainedTestFunction<Matrix, Vector> * fun_constrained = dynamic_cast<ConstrainedTestFunction<Matrix, Vector> *>(test_functions[i].get()); 
					solver.set_box_constraints(fun_constrained->box_constraints()); 
					solver.solve(*fun_constrained, x_init); 
					
					bool feas_flg = fun_constrained->is_feasible(x_init); 
					utopia_test_assert(feas_flg);

					// disp(x_init, "sol"); 
					// auto lb = fun_constrained->lower_bound(); 
					// auto ub = fun_constrained->upper_bound(); 
					// disp(lb, "lb"); 
					// disp(ub, "ub"); 
					// std::cout<<"--------------------------\n";
				}
				else
				{
					solver.solve(*test_functions[i], x_init); 
				}

				auto sol_status = solver.solution_status(); 
				//sol_status.describe(std::cout); 
				
				const auto dim = test_functions[i]->dim(); 
				const auto num_its = sol_status.iterates; 
				const auto conv_reason = sol_status.reason; 

				// utopia_test_assert(conv_reason > 0);

				if(exp_verbose && mpi_world_rank()==0)
				{
					std::cout<< i <<std::setw(5-std::to_string(i).size()) <<" : "<< test_functions[i]->name() <<"_" << dim <<  std::right <<  std::setw(60-std::to_string(dim).size() - test_functions[i]->name().size())  << std::right << "its:  " << num_its << std::setw(5-std::to_string(num_its).size())<<  "  \n"; 
				
					if(conv_reason< 0)
					{
						sol_status.describe(std::cout); 
					}

				}

				// if(test_functions[i]->exact_sol_known())
				// {
					// disp(x_init, "num_sol..."); 
					// disp(test_functions[i]->exact_sol(), "exact solution"); 

					// disp(test_functions[i]->exact_sol());
					// std::cout<<"norm(diff): "<< norm_infty(x_init - test_functions[i]->exact_sol()) << " \n"; 
					// disp(x_init);
					// disp(test_functions[i]->exact_sol());
				// }

			}
		}


	private:
		std::vector<std::shared_ptr<UnconstrainedExtendedTestFunction<Matrix, Vector> > >  test_functions_parallel_;
		std::vector<std::shared_ptr<UnconstrainedTestFunction<Matrix, Vector> > >  test_functions_;

		std::vector<std::shared_ptr<ConstrainedTestFunction<Matrix, Vector> > >  test_functions_constrained_;

		SizeType n_; 
		bool verbose_; 
	};

	static void quasi_tr()
	{
		int verbosity_level = 1;
		const int n_global = 10; 
		bool alg_verbose = true; 

		if(Utopia::instance().verbose()) {
			verbosity_level = 2;
		}

		#ifdef WITH_PETSC
			QuasiTRTest<PetscMatrix, PetscVector> bench1(n_global, alg_verbose);
			bench1.set_verbosity_level(verbosity_level);
			bench1.run();
		#endif //WITH_PETSC

	}

	UTOPIA_REGISTER_TEST_FUNCTION(quasi_tr);
}

