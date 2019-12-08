#include "utopia_Testing.hpp"
#include "utopia_Chrono.hpp"
#include "utopia_MPI.hpp"
#include "utopia.hpp"
#include "utopia_Base.hpp"
#include "utopia_Benchmark.hpp"
#include "utopia_LargeScaleIncludes.hpp"
#include "utopia_ConstrainedBenchmark.hpp"
#include <string>
#include <cassert>

namespace utopia 
{

	template<class Matrix, class Vector>
	class QPConstrainedBenchmark final: public Benchmark 
	{
	public:
		DEF_UTOPIA_SCALAR(Vector);

		std::string name() override
		{
			return "QPConstrainedBenchmark benchmark.";
		}

		QPConstrainedBenchmark(const SizeType & n = 10, const bool verbose = false): n_(n), verbose_(verbose)
		{
			// if(mpi_world_size()==1)
			// {
				// test_functions_parallel_.resize(2);
				// test_functions_parallel_[0] = std::make_shared<Poisson1D<Matrix, Vector> >(n_);
				// test_functions_parallel_[1] = std::make_shared<Bratu1D<Matrix, Vector> >(n_);

				// test_functions_.resize(3);
				// test_functions_[0] = std::make_shared<Rosenbrock01<Matrix, Vector> >();
				// test_functions_[1] = std::make_shared<QPTestFunction_2D<Matrix, Vector> >();
				// test_functions_[2] = std::make_shared<Watson20<Matrix, Vector> >();
			// }
			// else
			// {
				// test_functions_parallel_.resize(2);
				// test_functions_parallel_[0] = std::make_shared<Poisson1D<Matrix, Vector> >(n_*mpi_world_size());
				// test_functions_parallel_[1] = std::make_shared<Bratu1D<Matrix, Vector> >(n_*mpi_world_size());

				test_functions_.resize(1); 
				// test_functions_[0] = std::make_shared<Poisson1D<Matrix, Vector> >(n_, 4);
				test_functions_[0] = std::make_shared<Poisson2D<Matrix, Vector> >(n_);
			// }
		}

		~QPConstrainedBenchmark()
		{
			test_functions_.clear(); 
		}


		void initialize() override
		{

			this->register_experiment("MPGRP_Test",
				[this]() {
					MPGRP<Matrix, Vector> solver;

		            solver.verbose(true);
		            run_test(this->test_functions_, solver, "MPGRP_Test", this->verbose_);
				}
			);
			

		}

	private:

		template<class Fun, class NonlinearSolver>
		static void run_test(std::vector<std::shared_ptr<Fun> > & test_functions, NonlinearSolver &solver, const std::string & solv_name,  const bool & exp_verbose = false) 
		{

			InputParameters in;
			in.set("atol", 1e-6);
			in.set("rtol", 1e-11);
			in.set("stol", 1e-14);
			in.set("stol", 1e-14);
			in.set("delta_min", 1e-13); 
			in.set("max-it", 10000); 
			in.set("verbose", true);
			solver.read(in); 



			if(exp_verbose && mpi_world_rank()==0)
			{
				std::cout<<"--------------------------------------------------------- \n";
				std::cout<<"				" << solv_name << "				\n";
				std::cout<<"--------------------------------------------------------- \n";
			}

			std::cout<<"------------------ 0 \n"; 
			std::cout<<"test_functions.size() "<< test_functions.size() << " \n"; 


	    	for(size_t i =0; i < test_functions.size(); i++)
	    	{
				Vector x_init = test_functions[i]->initial_guess(); 
				// solver.solve(*test_functions[i], x_init); 

				std::cout<<"------------------ 1 \n"; 

				Vector g; 
				Matrix H; 
				test_functions[i]->gradient(x_init, g); 
				g *= -1.0; 
				test_functions[i]->hessian(x_init, H); 

				// disp(g); 
				// disp(H); 
				// exit(0);


				solver.set_box_constraints(test_functions[i]->box_constraints()); 
				solver.solve(H, g, x_init);

				std::cout<<"size_x: "<< size(x_init).get(0) << "  \n"; 

				// disp(x_init);


				Poisson2D<Matrix, Vector> * fun_poisson2D = dynamic_cast<Poisson2D<Matrix, Vector> *>(test_functions[i].get());
				fun_poisson2D->output_to_VTK(x_init, "Poisson2D_new.vtk");


				// auto sol_status = solver.solution_status(); 


				//sol_status.describe(std::cout); 
				
				// const auto dim = test_functions[i]->dim(); 
				// const auto num_its = sol_status.iterates; 
				// const auto conv_reason = sol_status.reason; 

				// utopia_test_assert(conv_reason > 0);

				// if(exp_verbose && mpi_world_rank()==0)
				// {
				// 	std::cout<< i <<std::setw(5-std::to_string(i).size()) <<" : "<< test_functions[i]->name() <<"_" << dim <<  std::right <<  std::setw(60-std::to_string(dim).size() - test_functions[i]->name().size())  << std::right << "its:  " << num_its << std::setw(5-std::to_string(num_its).size())<<  "  \n"; 
				
					// if(conv_reason< 0)
					// {
					// 	sol_status.describe(std::cout); 
					// }

				// }

				// if(test_functions[i]->exact_sol_known())
				// {
					// disp(x_init, "num_sol..."); 
					// disp(test_functions[i]->exact_sol(), "exact solution"); 
					// disp(x_init, "sol");
					// std::cout<<"norm(diff): "<< norm_infty(x_init - test_functions[i]->exact_sol()) << " \n"; 

				// }

			}
		}


	private:
		std::vector<std::shared_ptr<ConstrainedExtendedTestFunction<Matrix, Vector> > >    test_functions_;

		SizeType n_; 
		bool verbose_; 

	};

	static void qp_constrained()
	{
		int verbosity_level = 1;
		const int n_global = 50; 
		bool alg_verbose = true; 

		if(Utopia::instance().verbose()) {
			verbosity_level = 2;
		}

		#ifdef WITH_PETSC
			QPConstrainedBenchmark<PetscMatrix, PetscVector> bench1(n_global, alg_verbose);
			bench1.set_verbosity_level(verbosity_level);
			bench1.run();
		#endif //WITH_PETSC

	}

	UTOPIA_REGISTER_TEST_FUNCTION(qp_constrained);
}

