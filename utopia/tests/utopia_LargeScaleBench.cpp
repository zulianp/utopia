#include "utopia_Testing.hpp"
#include "utopia_Chrono.hpp"
#include "utopia_MPI.hpp"
#include "utopia.hpp"
#include "utopia_Benchmark.hpp"
#include "utopia_LargeScaleIncludes.hpp"
#include <string>
#include <cassert>

namespace utopia 
{

	template<class Matrix, class Vector>
	class LargeScaleUnconstrainedBenchmark : public Benchmark 
	{
	public:
		DEF_UTOPIA_SCALAR(Vector);

		virtual std::string name() override
		{
			return "LargeScaleUnconstrainedBenchmark benchmark.";
		}

		LargeScaleUnconstrainedBenchmark(const SizeType & n = 10, const bool verbose = false): n_(n), verbose_(verbose)
		{
			test_functions_.resize(1);
			// test_functions_[0] = std::make_shared<Bratu2D<Matrix, Vector> >(n_);

			// test_functions_[0] = std::make_shared<Poisson3D<Matrix, Vector> >(n_);
			// test_functions_[0] = std::make_shared<Morebv1D<Matrix, Vector> >(n_);

			// test_functions_[0] = std::make_shared<Poisson2D<Matrix, Vector> >(n_);

			// test_functions_[0] = std::make_shared<Poisson1D<Matrix, Vector> >(n_);
			// test_functions_[1] = std::make_shared<Bratu1D<Matrix, Vector> >(n_);


			test_functions_[0] = std::make_shared<NonEllipse2D<Matrix, Vector> >(n_);
			

			// auto fun = Poisson1D<Matrix, Vector>(n_); 
		}

		~LargeScaleUnconstrainedBenchmark()
		{
			test_functions_.clear(); 
		}


		void initialize() override
		{

			this->register_experiment("NewtonTest_FACTORIZATION",
				[this]() {
		            auto lin_solver = std::make_shared<utopia::Factorization<Matrix, Vector> >();
		            Newton<Matrix, Vector> solver(lin_solver);
		            solver.verbose(false);
		            run_tr(this->test_functions_, solver, "NewtonTest_FACTORIZATION", this->verbose_);
				}
			);

			// this->register_experiment("NewtonTest_CG_HOMEMADE_jacobi",
			// 	[this]() {
		 //            auto lin_solver = std::make_shared<utopia::ConjugateGradient<Matrix, Vector, HOMEMADE> >();
		 //            lin_solver->set_preconditioner(std::make_shared<GaussSeidel<Matrix, Vector> >());

		 //            lin_solver->atol(1e-10); 
		 //            lin_solver->max_it(500);
		 //            lin_solver->verbose(false); 

		 //            Newton<Matrix, Vector> solver(lin_solver);
		 //            run_tr(this->test_functions_, solver, "NewtonTest_CG_HOMEMADE_jacobi", this->verbose_);
			// 	}
			// );

			// this->register_experiment("NewtonTest_CG_PETSC_jacobi",
			// 	[this]() {
		 //            auto lin_solver = std::make_shared<utopia::ConjugateGradient<Matrix, Vector> >();
		 //            lin_solver->pc_type("sor"); 
   //                 // lin_solver->set_preconditioner(std::make_shared<GaussSeidel<Matrix, Vector> >());
		 //            // lin_solver->set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector> >());

		 //            lin_solver->atol(1e-10); 
		 //            lin_solver->max_it(500);
		 //            lin_solver->verbose(false); 

		 //            Newton<Matrix, Vector> solver(lin_solver);
		 //            run_tr(this->test_functions_, solver, "NewtonTest_CG_PETSC_inv_diag", this->verbose_);
			// 	}
			// );	

			// this->register_experiment("NewtonTest_STCG_inv_diag",
			// 	[this]() {
		 //            auto lin_solver = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();
		 //            lin_solver->set_preconditioner(std::make_shared<PointJacobi<Matrix, Vector> >());

		 //            lin_solver->atol(1e-10); 
		 //            lin_solver->max_it(300);

		 //            Newton<Matrix, Vector> solver(lin_solver);
		 //            run_tr(this->test_functions_, solver, "NewtonTest_STCG_inv_diag", this->verbose_);
			// 	}
			// );					

			// this->register_experiment("NewtonTest_BiCGStab_GS_PETSC",
			// 	[this]() {
		 //            auto lin_solver = std::make_shared<utopia::BiCGStab<Matrix, Vector> >();
		 //            lin_solver->set_preconditioner(std::make_shared<GaussSeidel<Matrix, Vector> >());

		 //            lin_solver->atol(1e-10); 
		 //            lin_solver->max_it(300);

		 //            Newton<Matrix, Vector> solver(lin_solver);
		 //            run_tr(this->test_functions_, solver, "NewtonTest_BiCGStab_GS", this->verbose_);
			// 	}
			// );		

			// this->register_experiment("NewtonTest_BiCGStab_HOMEMADE_GS",
			// 	[this]() {
		 //            auto lin_solver = std::make_shared<utopia::BiCGStab<Matrix, Vector, HOMEMADE> >();
		 //            lin_solver->set_preconditioner(std::make_shared<GaussSeidel<Matrix, Vector> >());

		 //            lin_solver->atol(1e-10); 
		 //            lin_solver->max_it(300);

		 //            Newton<Matrix, Vector> solver(lin_solver);
		 //            run_tr(this->test_functions_, solver, "NewtonTest_BiCGStab_HOMEMADE_GS", this->verbose_);
			// 	}
			// );		

			// #ifdef WITH_PETSC
			// 	this->register_experiment("NewtonTest_GMRES_backtracking",
			// 		[this]() {
			//             auto lin_solver = std::make_shared<utopia::GMRES<Matrix, Vector> >();
			//             lin_solver->set_preconditioner(std::make_shared<GaussSeidel<Matrix, Vector> >());

			//             lin_solver->atol(1e-10); 
			//             lin_solver->max_it(300);

			//             auto ls_strat  = std::make_shared<utopia::Backtracking<Vector> >();
			//             Newton<Matrix, Vector> solver(lin_solver);
			//             solver.set_line_search_strategy(ls_strat);		

			//             run_tr(this->test_functions_, solver, "NewtonTest_BiCGStab_HOMEMADE_backtracking", this->verbose_);
			// 		}
			// 	);		



			// 	this->register_experiment("NewtonTest_GMRES_SimpleBacktracking",
			// 		[this]() {
			//             auto lin_solver = std::make_shared<utopia::BiCGStab<Matrix, Vector> >();
			//             lin_solver->set_preconditioner(std::make_shared<GaussSeidel<Matrix, Vector> >());

			//             lin_solver->atol(1e-10); 
			//             lin_solver->max_it(300);

			//             auto ls_strat  = std::make_shared<utopia::SimpleBacktracking<Vector> >();            
			//             Newton<Matrix, Vector> solver(lin_solver);
			//             solver.set_line_search_strategy(ls_strat);	

			//             run_tr(this->test_functions_, solver, "NewtonTest_BiCGStab_HOMEMADE_SimpleBacktracking", this->verbose_);
			// 		}
			// 	);		
			// #endif


			// this->register_experiment("TR_STCG",
			// 	[this]() {
			// 		auto subproblem = std::make_shared<SteihaugToint<Matrix, Vector> >();
			// 		subproblem->pc_type("asm");
			// 		TrustRegion<Matrix, Vector> solver(subproblem);
			// 		solver.delta0(1e10);
			// 		run_tr(this->test_functions_, solver, "TR_STCG", this->verbose_);
			// 	}
			// );				

		}

	private:

		template<class NonlinearSolver>
		static void run_tr(std::vector<std::shared_ptr<UnconstrainedExtendedTestFunction<Matrix, Vector> > > & test_functions, NonlinearSolver &solver, const std::string & solv_name,  const bool & exp_verbose = false) 
		{

			InputParameters in;
			in.set("atol", 1e-7);
			in.set("rtol", 1e-11);
			in.set("stol", 1e-14);
			in.set("stol", 1e-14);
			in.set("delta_min", 1e-13); 
			in.set("max-it", 500); 
			in.set("verbose", false);
			solver.read(in); 



			if(exp_verbose && mpi_world_rank()==0)
			{
				std::cout<<"--------------------------------------------------------- \n";
				std::cout<<"				" << solv_name << "				\n";
				std::cout<<"--------------------------------------------------------- \n";
			}

	    	// for(size_t i =0; i < test_functions.size(); i++)
	    	for(auto i =0; i < 1; i++)
	    	{
				Vector x_init = test_functions[i]->initial_guess(); 
				solver.solve(*test_functions[i], x_init); 

				auto sol_status = solver.solution_status(); 
				//sol_status.describe(std::cout); 
				
				const auto dim = test_functions[i]->dim(); 
				const auto num_its = sol_status.iterates; 
				// const auto conv_reason = sol_status.reason; 

				if(exp_verbose && mpi_world_rank()==0)
				{
					std::cout<< i <<std::setw(5-std::to_string(i).size()) <<" : "<< test_functions[i]->name() <<"_" << dim <<  std::right <<  std::setw(60-std::to_string(dim).size() - test_functions[i]->name().size())  << std::right << "its:  " << num_its << std::setw(5-std::to_string(num_its).size())<<  "  \n"; 
				
					// if(conv_reason< 0)
					// {
					// 	sol_status.describe(std::cout); 
					// }
				}

				// disp(x_init);
				// Poisson3D<Matrix, Vector> * fun_bratu = dynamic_cast<Poisson3D<Matrix, Vector> *>(test_functions.back().get());
				// // fun_bratu->output_to_VTK(test_functions[i]->exact_sol(), "Poisson2D_exact.vtk");
				// fun_bratu->output_to_VTK(x_init, "Poisson3D.vtk");
				// fun_bratu->output_to_VTK(test_functions[i]->exact_sol(), "Poisson3D_exact.vtk");


				// if(test_functions[i]->exact_sol_known())
				// {
				// 	// disp(x_init, "num_sol..."); 
				// 	// disp(test_functions[i]->exact_sol(), "exact solution"); 
				// 	std::cout<<"norm(diff): "<< norm_infty(x_init - test_functions[i]->exact_sol()) << " \n"; 
				// 	// disp(x_init);
				// }

				mpi_world_barrier();
			}
		}


	private:
		std::vector<std::shared_ptr<UnconstrainedExtendedTestFunction<Matrix, Vector> > >  test_functions_;
		SizeType n_; 
		bool verbose_; 


	};

	static void unconstrained_large_scale()
	{
		int verbosity_level = 1;
		const int n_global = 10; 
		bool alg_verbose = false; 

		if(Utopia::instance().verbose()) {
			verbosity_level = 2;
		}

		#ifdef WITH_PETSC
			LargeScaleUnconstrainedBenchmark<PetscMatrix, PetscVector> bench1(n_global, alg_verbose);
			bench1.set_verbosity_level(verbosity_level);
			bench1.run();
		#endif //WITH_PETSC

	}

	UTOPIA_REGISTER_TEST_FUNCTION(unconstrained_large_scale);
}

