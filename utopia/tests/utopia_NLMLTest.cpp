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
	class RMTR_test : public Benchmark
	{
	public:
		DEF_UTOPIA_SCALAR(Vector);

		virtual std::string name() override
		{
			return "RMTR_test benchmark.";
		}

		RMTR_test(const SizeType & n = 6, const bool verbose = false): n_(n), n_levels_(3), verbose_(verbose)
		{
			// ml_problems_.resize(4);
			// // ml_problems_[0] =  std::make_shared<PetscMultilevelTestProblem<Matrix, Vector, Poisson2D<Matrix, Vector> > > (2, n_levels_, n_);
			// ml_problems_[0] =  std::make_shared<PetscMultilevelTestProblem<Matrix, Vector, Poisson3D<Matrix, Vector> > > (3, n_levels_, n_);
			// // ml_problems_[0] = std::make_shared<MultiLevelTestProblem1D<Matrix, Vector, Poisson1D<Matrix, Vector> > > (n_levels_, n_);

			// ml_problems_[1] =  std::make_shared<MultiLevelTestProblem1D<Matrix, Vector, Morebv1D<Matrix, Vector> > > (n_levels_, n_);

			// ml_problems_[2] =  std::make_shared<MultiLevelTestProblem1D<Matrix, Vector, Bratu1D<Matrix, Vector> > > (n_levels_, n_);
			// // ml_problems_[0] =  std::make_shared<PetscMultilevelTestProblem<Matrix, Vector, Bratu2D<Matrix, Vector> > > (2, n_levels_, n_);

			// ml_problems_[3] =  std::make_shared<PetscMultilevelTestProblem<Matrix, Vector, NonEllipse2D<Matrix, Vector> > > (2, n_levels_, n_);

			ml_problems_.resize(1);
			ml_problems_[0] = std::make_shared<MultiLevelTestProblem1D<Matrix, Vector, Poisson1D<Matrix, Vector> > > (n_levels_, n_);
			// ml_problems_[0] =  std::make_shared<PetscMultilevelTestProblem<Matrix, Vector, Poisson2D<Matrix, Vector> > > (2, n_levels_, n_);
		}

		~RMTR_test()
		{
			ml_problems_.clear();
		}


		void initialize() override
		{
			// this->register_experiment("RMTR_first_order_test",
			// 	[this]() {
		 //            // auto tr_strategy_coarse = std::make_shared<utopia::KSP_TR<Matrix, Vector> >("stcg", "lu", true);
		 //            // auto tr_strategy_fine = std::make_shared<utopia::Lanczos<Matrix, Vector> >("sor");

			// 		auto tr_strategy_coarse = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();
			// 		// tr_strategy_coarse->set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector> >());
			// 		tr_strategy_coarse->set_preconditioner(std::make_shared<IdentityPreconditioner<Vector> >());
			// 		tr_strategy_coarse->atol(1e-12);

		 //            auto tr_strategy_fine = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();
		 //            tr_strategy_fine->set_preconditioner(std::make_shared<IdentityPreconditioner<Vector> > ());
		 //            tr_strategy_fine->atol(1e-12);

		 //            auto rmtr = std::make_shared<RMTR_l2<Matrix, Vector, FIRST_ORDER> >(n_levels_);

		 //            // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
		 //            rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);
		 //            // rmtr->norm_schedule(MultilevelNormSchedule::OUTER_CYCLE);

		 //            // Set TR-QP strategies
		 //            rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
		 //            rmtr->set_fine_tr_strategy(tr_strategy_fine);

		 //            run_test(this->ml_problems_, rmtr, "RMTR_first_order_test", this->verbose_);
			// 	}
			// );

			// this->register_experiment("RMTR_first_order_MGOPT_test",
			// 	[this]() {
		 //            auto tr_strategy_coarse = std::make_shared<utopia::KSP_TR<Matrix, Vector> >("stcg", "lu", true);
		 //            // auto tr_strategy_fine = std::make_shared<utopia::Lanczos<Matrix, Vector> >("sor");

		 //            auto tr_strategy_fine = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();
		 //            tr_strategy_fine->atol(1e-12);


		 //            auto rmtr = std::make_shared<RMTR_l2<Matrix, Vector, FIRST_ORDER_MGOPT> >(n_levels_);

		 //            // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
		 //            rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);


		 //            // Set TR-QP strategies
		 //            rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
		 //            rmtr->set_fine_tr_strategy(tr_strategy_fine);

		 //            run_test(this->ml_problems_, rmtr, "RMTR_first_order_MGOPT_test", this->verbose_);
			// 	}
			// );

			// this->register_experiment("RMTR_second_order_test",
			// 	[this]() {
		 //            // auto tr_strategy_coarse = std::make_shared<utopia::KSP_TR<Matrix, Vector> >("stcg", "lu", true);
		 //            // // auto tr_strategy_fine = std::make_shared<utopia::Lanczos<Matrix, Vector> >("sor");
		 //            auto tr_strategy_fine = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();
		 //            // tr_strategy_fine->set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector> >());
		 //            tr_strategy_fine->set_preconditioner(std::make_shared<IdentityPreconditioner<Vector> >());
		 //            tr_strategy_fine->atol(1e-12);


			// 		auto tr_strategy_coarse = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();
			// 		// tr_strategy_coarse->set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector> >());
			// 		tr_strategy_coarse->set_preconditioner(std::make_shared<IdentityPreconditioner<Vector> >());
			// 		tr_strategy_coarse->atol(1e-12);


		 //            auto rmtr = std::make_shared<RMTR_l2<Matrix, Vector, SECOND_ORDER> >(n_levels_);

		 //            // Set TR-QP strategies
		 //            // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
		 //            rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);

		 //            rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
		 //            rmtr->set_fine_tr_strategy(tr_strategy_fine);

		 //            run_test(this->ml_problems_, rmtr, "RMTR_second_order_test", this->verbose_);
			// 	}
			// );

			// this->register_experiment("RMTR_galerkin_test",
			// 	[this]() {
		 //            // auto tr_strategy_coarse = std::make_shared<utopia::KSP_TR<Matrix, Vector> >("stcg", "lu", true);
		 //            // // auto tr_strategy_fine = std::make_shared<utopia::Lanczos<Matrix, Vector> >("sor");
		 //            // auto tr_strategy_fine = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();

		 //           	auto tr_strategy_fine = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();
		 //            // tr_strategy_fine->set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector> >());
		 //            tr_strategy_fine->set_preconditioner(std::make_shared<IdentityPreconditioner<Vector> >());
		 //            tr_strategy_fine->atol(1e-12);


			// 		auto tr_strategy_coarse = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();
			// 		// tr_strategy_coarse->set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector> >());
			// 		tr_strategy_coarse->set_preconditioner(std::make_shared<IdentityPreconditioner<Vector> >());
			// 		tr_strategy_coarse->atol(1e-12);


		 //            auto rmtr = std::make_shared<RMTR_l2<Matrix, Vector, GALERKIN> >(n_levels_);

		 //            // Set TR-QP strategies
		 //            // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
		 //            rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);
		 //            rmtr->norm_schedule(MultilevelNormSchedule::OUTER_CYCLE);

		 //            rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
		 //            rmtr->set_fine_tr_strategy(tr_strategy_fine);

		 //            run_test(this->ml_problems_, rmtr, "RMTR_galerkin_test", this->verbose_);
			// 	}
			// );




			this->register_experiment("RMTR_first_order_infty",
				[this]() {
		           	auto tr_strategy_fine = std::make_shared<utopia::MPGRP<Matrix, Vector> >();
		           	// auto tr_strategy_fine = std::make_shared<utopia::ProjectedGaussSeidel<Matrix, Vector> >();
		           	// tr_strategy_fine->use_symmetric_sweep(false);
		            tr_strategy_fine->atol(1e-12);
		            // tr_strategy_fine->verbose(true);

					auto tr_strategy_coarse = std::make_shared<utopia::MPGRP<Matrix, Vector> >();
					// auto tr_strategy_coarse = std::make_shared<utopia::ProjectedGaussSeidel<Matrix, Vector> >();
					tr_strategy_coarse->atol(1e-12);

		            // auto rmtr = std::make_shared<RMTR_inf<Matrix, Vector, TRBoundsGratton<Matrix, Vector>, FIRST_ORDER_MGOPT> >(n_levels_);
		            auto rmtr = std::make_shared<RMTR_inf<Matrix, Vector, TRGrattonBoxKornhuber<Matrix, Vector>, FIRST_ORDER> >(n_levels_);

		            // Set TR-QP strategies
		            // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
		            // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);

		            rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
		            rmtr->set_fine_tr_strategy(tr_strategy_fine);

		            run_test(this->ml_problems_, rmtr, "RMTR_first_order_infty", this->verbose_);
				}
			);

		}

	private:
		template<class Problem, class NonlinearSolver>
		static void run_test(std::vector<std::shared_ptr<Problem> > & ml_problems, NonlinearSolver &solver, const std::string & solv_name,  const bool & exp_verbose = false)
		{
			InputParameters in;
			in.set("atol", 1e-6);
			in.set("rtol", 1e-11);
			in.set("stol", 1e-14);
			in.set("stol", 1e-14);
			in.set("delta_min", 1e-13);
			in.set("max-it", 20);
			in.set("verbose", false);


            // RMTR specific parameters
            in.set("max_coarse_it", 2);
            in.set("max_sucessful_coarse_it", 1);
            in.set("max_QP_coarse_it", 1000);
            in.set("pre_smoothing_steps", 2);
            in.set("post_smoothing_steps", 2);
            in.set("max_sucessful_smoothing_it", 1);
            in.set("max_QP_smoothing_it", 8);
            // in.set("delta0", 0.001);
            in.set("delta0", 1e10);
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

	            // Transfers and objective functions
	            solver->set_transfer_operators(ml_problems[i]->get_transfer());
	            solver->set_functions(ml_problems[i]->get_functions());

	            // Solve
	            solver->solve(x);


				auto sol_status = solver->solution_status();

				// std::cout<<"it: "<< sol_status.iterates << "  \n";
				// std::cout<<"gradient_norm: "<< sol_status.gradient_norm << "  \n";
				// std::cout<<"reason: "<< sol_status.reason << "  \n";


				if(exp_verbose && mpi_world_rank()==0)
				{

		            if(UnconstrainedExtendedTestFunction<Matrix, Vector> * test_fun = dynamic_cast<UnconstrainedExtendedTestFunction<Matrix, Vector> *>(ml_problems[i]->get_functions().back().get()))
		            {
		            	const auto dim = test_fun->dim();
						// const auto num_its = sol_status.iterates;
						// const auto conv_reason = sol_status.reason;

						std::cout<< i <<std::setw(5-std::to_string(i).size()) <<" : "<< test_fun->name() << "   \n";

						// if(conv_reason< 0)
						// {
						// 	sol_status.describe(std::cout);
						// }
		            }
				}
			}
		}

	private:
		SizeType n_;
		SizeType n_levels_;
		bool verbose_;
		std::vector<std::shared_ptr<MultilevelTestProblemBase<Matrix, Vector> > > ml_problems_;
	};

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	template<class Matrix, class Vector>
	class QuasiRMTR_test : public Benchmark
	{
	public:
		DEF_UTOPIA_SCALAR(Vector);

		virtual std::string name() override
		{
			return "NLML_test benchmark.";
		}

		QuasiRMTR_test(const SizeType & n = 6, const bool verbose = false): n_(n), n_levels_(4), verbose_(verbose)
		{
			ml_problems_.resize(3);
			ml_problems_[0] = std::make_shared<MultiLevelTestProblem1D<Matrix, Vector, Poisson1D<Matrix, Vector> > > (n_levels_, n_);
			ml_problems_[1] =  std::make_shared<MultiLevelTestProblem1D<Matrix, Vector, Bratu1D<Matrix, Vector> > > (n_levels_, n_);
			ml_problems_[2] =  std::make_shared<PetscMultilevelTestProblem<Matrix, Vector, NonEllipse2D<Matrix, Vector> > > (2, n_levels_, n_);
		}

		~QuasiRMTR_test()
		{
			ml_problems_.clear();
		}


		void initialize() override
		{
			// this->register_experiment("RMTR_quasi_LBFGS_test",
			// 	[this]() {

		 //           	auto tr_strategy_fine = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();
		 //            tr_strategy_fine->set_preconditioner(std::make_shared<IdentityPreconditioner<Vector> >());
		 //            tr_strategy_fine->atol(1e-12);

			// 		auto tr_strategy_coarse = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();
			// 		tr_strategy_coarse->set_preconditioner(std::make_shared<IdentityPreconditioner<Vector> >());
			// 		tr_strategy_coarse->atol(1e-12);

		 //            auto rmtr = std::make_shared<QuasiRMTR<Matrix, Vector> >(n_levels_);

		 //            // Set TR-QP strategies
		 //            rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
		 //            // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);
		 //            // rmtr->norm_schedule(MultilevelNormSchedule::OUTER_CYCLE);

   //          		const SizeType memory_size = 5;
   //              	auto hess_approx   = std::make_shared<LBFGS<Vector> >(memory_size);
   //          		rmtr->set_hessian_approximation_strategy(hess_approx);


		 //            rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
		 //            rmtr->set_fine_tr_strategy(tr_strategy_fine);

		 //            run_test(this->ml_problems_, rmtr, "RMTR_quasi_LBFGS_test", this->verbose_);
			// 	}
			// );


			this->register_experiment("RMTR_quasi_LBFGS_test_inf_unconstrained",
				[this]() {

		           	auto tr_strategy_fine = std::make_shared<utopia::MPGRP<Matrix, Vector> >();
		            // tr_strategy_fine->set_preconditioner(std::make_shared<IdentityPreconditioner<Vector> >());
		            tr_strategy_fine->atol(1e-12);

					auto tr_strategy_coarse = std::make_shared<utopia::MPGRP<Matrix, Vector> >();
					// tr_strategy_coarse->set_preconditioner(std::make_shared<IdentityPreconditioner<Vector> >());
					tr_strategy_coarse->atol(1e-12);

		            auto rmtr = std::make_shared<QuasiRMTR_inf<Matrix, Vector, TRBoundsGratton<Matrix, Vector> > >(n_levels_);

		            // Set TR-QP strategies
		            rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
		            // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);
		            // rmtr->norm_schedule(MultilevelNormSchedule::OUTER_CYCLE);

            		const SizeType memory_size = 5;
                	auto hess_approx   = std::make_shared<LBFGS<Vector> >(memory_size);
            		rmtr->set_hessian_approximation_strategy(hess_approx);


		            rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
		            rmtr->set_fine_tr_strategy(tr_strategy_fine);

		            run_test(this->ml_problems_, rmtr, "RMTR_quasi_LBFGS_test", this->verbose_);
				}
			);


			// this->register_experiment("RMTR_quasi_LSR1_test",
			// 	[this]() {

		 //           	auto tr_strategy_fine = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();
		 //            tr_strategy_fine->set_preconditioner(std::make_shared<IdentityPreconditioner<Vector> >());
		 //            tr_strategy_fine->atol(1e-12);

			// 		auto tr_strategy_coarse = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();
			// 		tr_strategy_coarse->set_preconditioner(std::make_shared<IdentityPreconditioner<Vector> >());
			// 		tr_strategy_coarse->atol(1e-12);

		 //            auto rmtr = std::make_shared<QuasiRMTR<Matrix, Vector> >(n_levels_);

		 //            // Set TR-QP strategies
		 //            // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
		 //            rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);

   //              	auto hess_approx   = std::make_shared<LSR1<Vector> >(5);
   //          		rmtr->set_hessian_approximation_strategy(hess_approx);


		 //            rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
		 //            rmtr->set_fine_tr_strategy(tr_strategy_fine);

		 //            run_test(this->ml_problems_, rmtr, "RMTR_quasi_LSR1_test", this->verbose_);
			// 	}
			// );


		}

	private:
		template<class Problem, class NonlinearSolver>
		static void run_test(std::vector<std::shared_ptr<Problem> > & ml_problems, NonlinearSolver &solver, const std::string & solv_name,  const bool & exp_verbose = false)
		{
			InputParameters in;
			in.set("atol", 1e-6);
			in.set("rtol", 1e-11);
			in.set("stol", 1e-14);
			in.set("stol", 1e-14);
			in.set("delta_min", 1e-13);
			in.set("max-it", 50);
			in.set("verbose", false);

            // RMTR specific parameters
            in.set("max_coarse_it", 10);
            in.set("max_sucessful_coarse_it", 5);
            in.set("max_QP_coarse_it", 1000);
            in.set("pre_smoothing_steps", 10);
            in.set("post_smoothing_steps", 10);
            in.set("max_sucessful_smoothing_it", 5);
            // in.set("max_QP_smoothing_it", 10);
            in.set("max_QP_smoothing_it", 1000);
            // in.set("delta0", 1.0e10);
            in.set("grad_smoothess_termination", 1e-8);
            // in.set("skip_BC_checks", true);

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
				assert(!empty(x));

				// x.set(20);

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



	static void rmtr()
	{
		int verbosity_level = 0;
		const int n_global = 10;
		bool alg_verbose = false;

		//FIXME create a special purpose flag only for these
		// if(Utopia::instance().verbose()) {
		// 	verbosity_level = 2;
		// }

		#ifdef WITH_PETSC
			RMTR_test<PetscMatrix, PetscVector> bench1(n_global, alg_verbose);
			bench1.set_verbosity_level(verbosity_level);
			bench1.run();
		#endif //WITH_PETSC

	}

	static void quasi_rmtr()
	{
		int verbosity_level = 0;
		const int n_global = 10;
		bool alg_verbose = false;

		//FIXME create a special purpose flag only for these
		// if(Utopia::instance().verbose()) {
		// 	verbosity_level = 2;
		// }

		#ifdef WITH_PETSC
			QuasiRMTR_test<PetscMatrix, PetscVector> bench1(n_global, alg_verbose);
			bench1.set_verbosity_level(verbosity_level);
			bench1.run();
		#endif //WITH_PETSC

	}



	UTOPIA_REGISTER_TEST_FUNCTION(rmtr);
	UTOPIA_REGISTER_TEST_FUNCTION(quasi_rmtr);
}

