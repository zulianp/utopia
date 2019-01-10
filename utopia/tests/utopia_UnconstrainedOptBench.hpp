#ifndef UTOPIA_UNCONSTRAINED_OPTIMIZATION_BENCHMARK_HPP
#define UTOPIA_UNCONSTRAINED_OPTIMIZATION_BENCHMARK_HPP

#include "utopia_Chrono.hpp"
#include "utopia_MPI.hpp"
#include "utopia.hpp"
#include "utopia_Benchmark.hpp"
#include "utopia_assemble_laplacian_1D.hpp"


#include <string>
#include <cassert>

namespace utopia 
{

	template<class Matrix, class Vector>
	class UnconstrainedOptimizationBenchmark : public Benchmark 
	{
	public:
		DEF_UTOPIA_SCALAR(Vector)

		virtual std::string name() override
		{
			return "TR: unconstrained optimization benchmark";
		}

		UnconstrainedOptimizationBenchmark(const SizeType & n = 10) : n_(n), verbose_(true)
		{
			test_functions_.resize(18); 
			test_functions_[0] = std::make_shared<Hellical07<Matrix, Vector> >();
	    	test_functions_[1] = std::make_shared<Biggs18<Matrix, Vector> >();	
			test_functions_[2] = std::make_shared<Gaussian09<Matrix, Vector> >();	    			
	    	test_functions_[3] = std::make_shared<Powell03<Matrix, Vector> >(); 
	    	test_functions_[4] = std::make_shared<Box12<Matrix, Vector> >();	  
			
			test_functions_[5] = std::make_shared<VariablyDim25<Matrix, Vector> >(n_); 	 // works also in parallel 		
	    	test_functions_[6] = std::make_shared<Watson20<Matrix, Vector> >(); 		
			test_functions_[7] = std::make_shared<PenaltyI23<Matrix, Vector> >(n_); 	    // works also in parallel 		    			    	  	
	    	test_functions_[8] = std::make_shared<Brown04<Matrix, Vector> >();
	    	test_functions_[9] = std::make_shared<PenaltyII24<Matrix, Vector> >();

	    	test_functions_[10] = std::make_shared<BrownDennis16<Matrix, Vector> >();	
	    	test_functions_[11] = std::make_shared<Gulf11<Matrix, Vector> >(); 
			test_functions_[12] = std::make_shared<Trigonometric26<Matrix, Vector> >(n_);  	
	    	test_functions_[13] = std::make_shared<ExtendedRosenbrock21<Matrix, Vector> >(n_); // works also in parallel 		    			    	  	
	    	test_functions_[14] = std::make_shared<Beale05<Matrix, Vector> >();

	    	test_functions_[15] = std::make_shared<Woods14<Matrix, Vector> >();
	    	test_functions_[16] = std::make_shared<ExtendedPowell22<Matrix, Vector> >(12);
	    	test_functions_[17] = std::make_shared<Chebyquad35<Matrix, Vector> >();			
		}

		~UnconstrainedOptimizationBenchmark()
		{
			test_functions_.clear(); 
		}


		void initialize() override
		{

			this->register_experiment("TR_STCG",
				[this]() {
					auto subproblem = std::make_shared<SteihaugToint<Matrix, Vector> >();
					TrustRegion<Matrix, Vector> solver(subproblem);
					run_tr(this->test_functions_, solver, "TR_STCG", this->verbose_);
				}
			);

			this->register_experiment("TR_Lanczos",
				[this]() {
					auto subproblem = std::make_shared<Lanczos<Matrix, Vector> >();
					TrustRegion<Matrix, Vector> solver(subproblem);
					run_tr(this->test_functions_, solver, "TR_Lanczos", this->verbose_);
				}
			);			

			this->register_experiment("TR_Nash",
				[this]() {
					auto subproblem = std::make_shared<Nash<Matrix, Vector> >();
					TrustRegion<Matrix, Vector> solver(subproblem);
					run_tr(this->test_functions_, solver, "TR_Nash", this->verbose_);
				}
			);		

			this->register_experiment("TR_Dogleg",
				[this]() {
					auto linear_solver = std::make_shared<GMRES<Matrix, Vector>>();	
					linear_solver->atol(1e-14); 
					linear_solver->max_it(10000);
					auto subproblem = std::make_shared<Dogleg<Matrix, Vector> >(linear_solver); 
					TrustRegion<Matrix, Vector> solver(subproblem);
					run_tr(this->test_functions_, solver, "TR_Dogleg", this->verbose_);
				}
			);			

			// // TODO:: add check for slepcs
			// this->register_experiment("TR_MS",
			// 	[this]() {
			// 		auto eigen_solver = std::make_shared<SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL> >();
			// 		// TODO:: add checks if has arpack
			// 		eigen_solver->solver_type("arpack");
					
			// 		auto linear_solver = std::make_shared<LUDecomposition<Matrix, Vector> >();
			// 		linear_solver->set_library_type(PETSC_TAG); 

			// 		auto subproblem = std::make_shared<utopia::MoreSorensenEigen<DMatrixd, Vector> >(linear_solver, eigen_solver);
			// 		TrustRegion<Matrix, Vector> solver(subproblem);
			// 		run_tr(this->test_functions_, solver, "TR_MS", this->verbose_);
			// 	}
			// );						


		}

	private:
		SizeType n_; 
		std::vector<std::shared_ptr<UnconstrainedTestFunction<Matrix, Vector> > >  test_functions_;
		bool verbose_; 

		template<class TRSolver>
		static void run_tr(std::vector<std::shared_ptr<UnconstrainedTestFunction<Matrix, Vector> > > & test_functions, TRSolver &solver, const std::string & solv_name,  const bool & exp_verbose = false) 
		{

			InputParameters in;
			in.set("atol", 1e-7);
			in.set("rtol", 1e-11);
			in.set("stol", 1e-14);
			in.set("stol", 1e-14);
			in.set("delta_min", 1e-13); 
			in.set("max_it", 300); 
			in.set("verbose", false); 

			auto params_qp = std::make_shared<InputParameters>(); 
			params_qp->set("atol", 1e-14); 
			params_qp->set("rtol", 1e-14); 
			params_qp->set("stol", 1e-14); 
			auto params_qp_cast = std::static_pointer_cast<Input>(params_qp); 				
			in.set("linear-solver", params_qp_cast);


			solver.read(in); 

			if(exp_verbose){
				std::cout<<"--------------------------------------------------------- \n";
				std::cout<<"				" << solv_name << "				\n";
				std::cout<<"--------------------------------------------------------- \n";
			}

	    	for(auto i =0; i < test_functions.size(); i++)
	    	{
				DVectord x_init = test_functions[i]->initial_guess(); 
				solver.solve(*test_functions[i], x_init); 

				auto sol_status = solver.solution_status(); 
				//sol_status.describe(std::cout); 
				
				const auto dim = test_functions[i]->dim(); 
				const auto num_its = sol_status.iterates; 

				if(exp_verbose)
				{
					std::cout<< i <<std::setw(5-std::to_string(i).size()) <<" : "<< test_functions[i]->name() <<"_" << dim <<  std::right <<  std::setw(40-std::to_string(dim).size() - test_functions[i]->name().size())  << std::right << "its:  " << num_its << "  \n"; 
				}

				if(test_functions[i]->exact_sol_known())
				{
					// disp(x_init); 
					utopia_test_assert(approxeq(x_init, test_functions[i]->exact_sol(), 1e-4));
				}
			}
		}
	};


	void run_unconstrained_optimization_benchmark()
	{
		int verbosity_level = 1;
		if(Utopia::instance().verbose()) {
			verbosity_level = 2;
		}

		if(mpi_world_size()==1)
		{
			#ifdef WITH_PETSC
				UnconstrainedOptimizationBenchmark<DMatrixd, DVectord> bench1;
				bench1.set_verbosity_level(verbosity_level);
				bench1.run();
			#endif //WITH_PETSC
		}
		else
		{
			std::cout<<"run_unconstrained_optimization_benchmark, does not work in parallel. \n"; 
		}
	}


}

#endif //UTOPIA_UNCONSTRAINED_OPTIMIZATION_BENCHMARK_HPP
