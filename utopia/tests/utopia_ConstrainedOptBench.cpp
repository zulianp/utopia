#include "utopia_Testing.hpp"
#include "utopia_Chrono.hpp"
#include "utopia_MPI.hpp"
#include "utopia.hpp"
#include "utopia_Benchmark.hpp"
#include "utopia_assemble_laplacian_1D.hpp"
#include "utopia_ConstrainedBenchmark.hpp"
#include <string>
#include <cassert>

namespace utopia 
{

	template<class Matrix, class Vector>
	class ConstrainedOptimizationBenchmark : public Benchmark 
	{
	public:
		DEF_UTOPIA_SCALAR(Vector);

		virtual std::string name() override
		{
			return "TR: Bound constrained optimization benchmark";
		}

		ConstrainedOptimizationBenchmark(): verbose_(true)
		{
			test_functions_.resize(18); 
			test_functions_[0] = std::make_shared<Powell03Constrained<Matrix, Vector> >();
			test_functions_[1] = std::make_shared<Beale05Constrained<Matrix, Vector> >();
			test_functions_[2] = std::make_shared<Hellical07Constrained<Matrix, Vector> >();
			test_functions_[3] = std::make_shared<Gaussian09Constrained<Matrix, Vector> >();
			test_functions_[4] = std::make_shared<Gulf11Constrained<Matrix, Vector> >();

			test_functions_[5] = std::make_shared<Box12Constrained<Matrix, Vector> >();
			test_functions_[6] = std::make_shared<Woods14Constrained<Matrix, Vector> >();
			test_functions_[7] = std::make_shared<BrownDennis16Constrained<Matrix, Vector> >();
			test_functions_[8] = std::make_shared<Biggs18Constrained<Matrix, Vector> >();
			test_functions_[9] = std::make_shared<Rosenbrock21Constrained<Matrix, Vector> >();
			
			test_functions_[10] = std::make_shared<ExtendedPowell22Constrained<Matrix, Vector> >();
			test_functions_[11] = std::make_shared<PenaltyI23Constrained<Matrix, Vector> >();
			test_functions_[12] = std::make_shared<VariablyDim25Constrained<Matrix, Vector> >();
			test_functions_[13] = std::make_shared<Trigonometric26Constrained<Matrix, Vector> >();
			test_functions_[14] = std::make_shared<Chebyquad35Constrained<Matrix, Vector> >(); 
			
			test_functions_[15] = std::make_shared<Brown04Constrained<Matrix, Vector> >(); 
			test_functions_[16] = std::make_shared<Watson20Constrained<Matrix, Vector> >(); 
			test_functions_[17] = std::make_shared<PenaltyII24Constrained<Matrix, Vector> >();

		}

		~ConstrainedOptimizationBenchmark()
		{
			test_functions_.clear(); 
		}


		void initialize() override
		{

			this->register_experiment("TR_Variable_Bound_MPRGP",
				[this]() {
		            auto subproblem = std::make_shared<utopia::MPGRP<Matrix, Vector> >();
		            subproblem->atol(1e-14); 
		            subproblem->stol(1e-14); 
		            subproblem->rtol(1e-14); 
		            subproblem->verbose(false);

		            TrustRegionVariableBound<Matrix, Vector> tr_solver(subproblem);
		            run_tr(this->test_functions_, tr_solver, "TR_Variable_Bound_MPRGP", this->verbose_);
				}
			);

			this->register_experiment("TR_Variable_BGS",
				[this]() {
		            auto subproblem = std::make_shared<utopia::ProjectedGaussSeidel<Matrix, Vector> >();
		            subproblem->atol(1e-14); 
		            subproblem->stol(1e-14); 
		            subproblem->rtol(1e-14); 
		            subproblem->verbose(false);
		            subproblem->use_line_search(false); 

		            TrustRegionVariableBound<Matrix, Vector> tr_solver(subproblem);
		            run_tr(this->test_functions_, tr_solver, "TR_Variable_BGS", this->verbose_);
				}
			);	

		#ifdef WITH_PETSC
			this->register_experiment("TR_Variable_Tao",
				[this]() {
		            auto lsolver = std::make_shared<LUDecomposition<PetscMatrix, PetscVector> >();
		            lsolver->set_library_type("petsc"); 
		            auto subproblem =  std::make_shared<utopia::TaoQPSolver<PetscMatrix, PetscVector> >(lsolver);

		            TrustRegionVariableBound<Matrix, Vector> tr_solver(subproblem);
		            run_tr(this->test_functions_, tr_solver, "TR_Variable_Tao", this->verbose_);
				}
			);		
		#endif //WITH_PETSC							


			this->register_experiment("TR_Variable_ProjGradient",
				[this]() {
		            auto subproblem = std::make_shared<utopia::ProjectedGradient<Matrix, Vector> >();
		            subproblem->atol(1e-14); 
		            subproblem->stol(1e-14); 
		            subproblem->rtol(1e-14); 
		            subproblem->verbose(false);

		            TrustRegionVariableBound<Matrix, Vector> tr_solver(subproblem);
		            run_tr(this->test_functions_, tr_solver, "TR_Variable_ProjGradient", this->verbose_);
				}
			);	


			this->register_experiment("TR_Variable_ProjCG",
				[this]() {
		            auto subproblem = std::make_shared<utopia::ProjectedConjugateGradient<Matrix, Vector> >();
		            subproblem->atol(1e-14); 
		            subproblem->stol(1e-14); 
		            subproblem->rtol(1e-14); 
		            subproblem->verbose(false);

		            TrustRegionVariableBound<Matrix, Vector> tr_solver(subproblem);
		            run_tr(this->test_functions_, tr_solver, "TR_Variable_ProjCG", this->verbose_);
				}
			);	


			this->register_experiment("TR_Variable_SemiSmoothNewton",
				[this]() {
					auto lsolver = std::make_shared<LUDecomposition<PetscMatrix, PetscVector> >();
		            auto subproblem = std::make_shared<utopia::SemismoothNewton<Matrix, Vector> >(lsolver);
		            subproblem->atol(1e-14); 
		            subproblem->stol(1e-14); 
		            subproblem->rtol(1e-14); 
		            subproblem->verbose(false);

		            TrustRegionVariableBound<Matrix, Vector> tr_solver(subproblem);
		            run_tr(this->test_functions_, tr_solver, "TR_Variable_SemiSmoothNewton", this->verbose_);
				}
			);	



		}

	private:
		std::vector<std::shared_ptr<ConstrainedTestFunction<Matrix, Vector> > >  test_functions_;
		bool verbose_; 

		template<class TRSolver>
		static void run_tr(std::vector<std::shared_ptr<ConstrainedTestFunction<Matrix, Vector> > > & test_functions, TRSolver &solver, const std::string & solv_name,  const bool & exp_verbose = false) 
		{

			InputParameters in;
			in.set("atol", 1e-6);
			in.set("rtol", 1e-11);
			in.set("stol", 1e-14);
			in.set("stol", 1e-14);
			in.set("delta_min", 1e-13); 
			in.set("max-it", 5000); 
			in.set("verbose", false); 

			// auto params_qp = std::make_shared<InputParameters>(); 
			// params_qp->set("atol", 1e-14); 
			// params_qp->set("rtol", 1e-14); 
			// params_qp->set("stol", 1e-14); 
			// auto params_qp_cast = std::static_pointer_cast<Input>(params_qp); 				
			// in.set("linear-solver", params_qp_cast);
			solver.read(in); 



			if(exp_verbose){
				std::cout<<"--------------------------------------------------------- \n";
				std::cout<<"				" << solv_name << "				\n";
				std::cout<<"--------------------------------------------------------- \n";
			}

	    	for(size_t i =0; i < test_functions.size(); i++)
	    	// for(auto i =0; i < 1; i++)
	    	{
				Vector x_init = test_functions[i]->initial_guess(); 

            	solver.set_box_constraints(test_functions[i]->box_constraints());
				solver.solve(*test_functions[i], x_init); 

				auto sol_status = solver.solution_status(); 
				//sol_status.describe(std::cout); 
				
				const auto dim = test_functions[i]->dim(); 
				const auto num_its = sol_status.iterates; 
				// const auto conv_reason = sol_status.reason; 

				if(exp_verbose)
				{
					std::cout<< i <<std::setw(5-std::to_string(i).size()) <<" : "<< test_functions[i]->name() <<"_" << dim <<  std::right <<  std::setw(60-std::to_string(dim).size() - test_functions[i]->name().size())  << std::right << "its:  " << num_its << std::setw(5-std::to_string(num_its).size())<<  "  \n"; 
				
					// if(conv_reason< 0)
					// {
					// 	sol_status.describe(std::cout); 
					// }
				}

				if(test_functions[i]->exact_sol_known())
				{
					// disp(x_init); 
				}
			}
		}
	};

	static void constrained_opt()
	{
		int verbosity_level = 1;
		if(Utopia::instance().verbose()) {
			verbosity_level = 2;
		}

		if(mpi_world_size()==1)
		{
			#ifdef WITH_PETSC
				ConstrainedOptimizationBenchmark<PetscMatrix, PetscVector> bench1;
				bench1.set_verbosity_level(verbosity_level);
				bench1.run();
			#endif //WITH_PETSC
		}
		else
		{
			std::cout<<"constrained_opt, does not work in parallel. \n"; 
		}
	}

	UTOPIA_REGISTER_TEST_FUNCTION_OPTIONAL(constrained_opt);
}

