#include "utopia_Testing.hpp"
#include "utopia_Chrono.hpp"
#include "utopia_MPI.hpp"
#include "utopia.hpp"
#include "utopia_Base.hpp"
#include "utopia_Benchmark.hpp"
#include "utopia_LargeScaleIncludes.hpp"
#include <string>
#include <cassert>

namespace utopia 
{

	template<class Matrix, class Vector>
	class NLML_test : public Benchmark 
	{
	public:
		DEF_UTOPIA_SCALAR(Vector);

		virtual std::string name() override
		{
			return "NLML_test benchmark.";
		}

		NLML_test(const SizeType & n = 5, const bool verbose = false): n_(n), n_levels_(5), verbose_(verbose)
		{
			ml_problems_.resize(1); 
			ml_problems_[0] =  std::make_shared<PetscMultilevelTestProblem<Matrix, Vector, Poisson2D<Matrix, Vector> > > (2, n_levels_, n_);

		}

		~NLML_test()
		{
			ml_problems_.clear(); 
		}


		void initialize() override
		{
			this->register_experiment("RMTR_first_order_test",
				[this]() {
		            auto tr_strategy_coarse = std::make_shared<utopia::KSP_TR<Matrix, Vector> >("stcg", "lu", true);
		            auto tr_strategy_fine = std::make_shared<utopia::Lanczos<Matrix, Vector> >("sor");
		            auto rmtr = std::make_shared<RMTR<Matrix, Vector, FIRST_ORDER> >(n_levels_);

		            // Set TR-QP strategies
		            rmtr->verbose(true);
		            rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
		            rmtr->norm_schedule(MultilevelNormSchedule::OUTER_CYCLE);

		            rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
		            rmtr->set_fine_tr_strategy(tr_strategy_fine);

		            run_test(this->ml_problems_, rmtr, "RMTR_first_order_test", this->verbose_);
				}
			);

		}

	private:
		template<class Problem, class NonlinearSolver>
		static void run_test(std::vector<std::shared_ptr<Problem> > & ml_problems, NonlinearSolver &solver, const std::string & solv_name,  const bool & exp_verbose = false) 
		{
			InputParameters in;
			in.set("atol", 1e-7);
			in.set("rtol", 1e-11);
			in.set("stol", 1e-14);
			in.set("stol", 1e-14);
			in.set("delta_min", 1e-13); 
			in.set("max-it", 50); 
			in.set("verbose", true);

            // RMTR specific parameters
            in.set("max_coarse_it", 1);
            in.set("max_sucessful_coarse_it", 1);
            in.set("max_QP_coarse_it", 1000);
            in.set("pre_smoothing_steps", 1);
            in.set("post_smoothing_steps", 1);
            in.set("max_sucessful_smoothing_it", 1);
            in.set("max_QP_smoothing_it", 1);
            in.set("delta0", 1.0e3);
            in.set("grad_smoothess_termination", 1e-8);

			solver->read(in); 

			if(exp_verbose && mpi_world_rank()==0)
			{
				std::cout<<"--------------------------------------------------------- \n";
				std::cout<<"				" << solv_name << "				\n";
				std::cout<<"--------------------------------------------------------- \n";
			}

	    	for(size_t i =0; i < ml_problems.size(); i++)
	    	{
				Vector x = ml_problems[i]->get_functions().back()->initial_guess(); 
				 
				x.set(20);

	            // Transfers and objective functions
	            solver->set_transfer_operators(ml_problems[i]->get_transfer());
	            solver->set_functions(ml_problems[i]->get_functions());

	            // Solve
	            solver->solve(x);


				// auto sol_status = solver.solution_status();  
				// const auto dim = test_functions[i]->dim(); 
				// const auto num_its = sol_status.iterates; 
				// const auto conv_reason = sol_status.reason; 

				// if(exp_verbose && mpi_world_rank()==0)
				// {
				// 	std::cout<< i <<std::setw(5-std::to_string(i).size()) <<" : "<< test_functions[i]->name() <<"_" << dim <<  std::right <<  std::setw(60-std::to_string(dim).size() - test_functions[i]->name().size())  << std::right << "its:  " << num_its << std::setw(5-std::to_string(num_its).size())<<  "  \n"; 
				
				// 	if(conv_reason< 0)
				// 	{
				// 		sol_status.describe(std::cout); 
				// 	}

				// }
			}
		}

	private:
		SizeType n_; 
		SizeType n_levels_; 
		bool verbose_; 
		std::vector<std::shared_ptr<MultilevelTestProblemBase<Matrix, Vector> > > ml_problems_; 
	};

	static void nlml_test()
	{
		int verbosity_level = 1;
		const int n_global = 10; 
		bool alg_verbose = true; 

		if(Utopia::instance().verbose()) {
			verbosity_level = 2;
		}

		#ifdef WITH_PETSC
			NLML_test<PetscMatrix, PetscVector> bench1(n_global, alg_verbose);
			bench1.set_verbosity_level(verbosity_level);
			bench1.run();
		#endif //WITH_PETSC

	}

	UTOPIA_REGISTER_TEST_FUNCTION(nlml_test);
}

