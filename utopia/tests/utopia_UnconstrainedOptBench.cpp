#include "utopia.hpp"
#include "utopia_Benchmark.hpp"
#include "utopia_Chrono.hpp"
#include "utopia_MPI.hpp"
#include "utopia_Testing.hpp"
#include "utopia_UnconstrainedBenchmark.hpp"
#include "utopia_assemble_laplacian_1D.hpp"

#ifdef WITH_PETSC
#ifdef WITH_SLEPC
#include "utopia_petsc_Slepc.hpp"
#endif  // WITH_SLEPC
#endif  // WITH_PETSC

#include <cassert>
#include <string>

namespace utopia {

    template <class Matrix, class Vector>
    class UnconstrainedOptimizationBenchmark : public Benchmark {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        std::string name() override { return "TR: unconstrained optimization benchmark"; }

        explicit UnconstrainedOptimizationBenchmark(const SizeType &n) : n_(n) {
            test_functions_.resize(18);
            test_functions_[0] = std::make_shared<Hellical07<Matrix, Vector>>();
            test_functions_[1] = std::make_shared<Biggs18<Matrix, Vector>>();
            test_functions_[2] = std::make_shared<Gaussian09<Matrix, Vector>>();
            test_functions_[3] = std::make_shared<Powell03<Matrix, Vector>>();
            test_functions_[4] = std::make_shared<Box12<Matrix, Vector>>();

            test_functions_[5] = std::make_shared<VariablyDim25<Matrix, Vector>>(n_);  // works also in parallel
            test_functions_[6] = std::make_shared<Watson20<Matrix, Vector>>();
            test_functions_[7] = std::make_shared<PenaltyI23<Matrix, Vector>>(n_);  // works also in parallel
            test_functions_[8] = std::make_shared<PenaltyII24<Matrix, Vector>>();
            test_functions_[9] = std::make_shared<Brown04<Matrix, Vector>>();

            test_functions_[10] = std::make_shared<BrownDennis16<Matrix, Vector>>();
            test_functions_[11] = std::make_shared<Gulf11<Matrix, Vector>>();
            test_functions_[12] = std::make_shared<Trigonometric26<Matrix, Vector>>(n_);
            test_functions_[13] = std::make_shared<ExtendedRosenbrock21<Matrix, Vector>>(n_);  // works also in parallel
            test_functions_[14] = std::make_shared<Beale05<Matrix, Vector>>();

            test_functions_[15] = std::make_shared<Woods14<Matrix, Vector>>();
            test_functions_[16] = std::make_shared<ExtendedPowell22<Matrix, Vector>>(64);
            test_functions_[17] = std::make_shared<Chebyquad35<Matrix, Vector>>();
        }

        ~UnconstrainedOptimizationBenchmark() override { test_functions_.clear(); }

        void initialize() override {
            this->register_experiment("TR_STCG_HOMEMADE", [this]() {
                auto subproblem = std::make_shared<SteihaugToint<Matrix, Vector, HOMEMADE>>();
                subproblem->set_preconditioner(std::make_shared<IdentityPreconditioner<Vector>>());
                TrustRegion<Matrix, Vector> solver(subproblem);
                run_tr(this->test_functions_, solver, "TR_STCG_HOMEMADE", this->verbose_);
            });

            this->register_experiment("TR_STCG", [this]() {
                auto subproblem = std::make_shared<SteihaugToint<Matrix, Vector>>();
                // subproblem->pc_type("none");
                TrustRegion<Matrix, Vector> solver(subproblem);
                run_tr(this->test_functions_, solver, "TR_STCG", this->verbose_);
            });
#ifdef WITH_PETSC
            this->register_experiment("TR_Lanczos", [this]() {
                auto subproblem = std::make_shared<Lanczos<Matrix, Vector>>();
                TrustRegion<Matrix, Vector> solver(subproblem);
                run_tr(this->test_functions_, solver, "TR_Lanczos", this->verbose_);
            });

            this->register_experiment("TR_Nash", [this]() {
                auto subproblem = std::make_shared<Nash<Matrix, Vector>>();
                TrustRegion<Matrix, Vector> solver(subproblem);
                run_tr(this->test_functions_, solver, "TR_Nash", this->verbose_);
            });

#endif  // WITH_PETSC

            this->register_experiment("TR_Dogleg", [this]() {
                auto linear_solver = std::make_shared<GMRES<Matrix, Vector>>();
                linear_solver->atol(1e-14);
                linear_solver->max_it(10000);
                auto subproblem = std::make_shared<Dogleg<Matrix, Vector>>(linear_solver);
                TrustRegion<Matrix, Vector> solver(subproblem);
                run_tr(this->test_functions_, solver, "TR_Dogleg", this->verbose_);
            });

            this->register_experiment("TR_Cauchy_Point", [this]() {
                auto subproblem = std::make_shared<CauchyPoint<Matrix, Vector>>();
                TrustRegion<Matrix, Vector> solver(subproblem);
                run_tr(this->test_functions_, solver, "TR_Cauchy_Point", this->verbose_);
            });

#ifdef WITH_PETSC
#ifdef WITH_SLEPC
            // TODO(zulianp): : add check for slepcs
            this->register_experiment("TR_MS", [this]() {
                auto eigen_solver = std::make_shared<SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL>>();
                // TODO(zulianp): : add checks if has arpack
                eigen_solver->solver_type("arpack");

                auto linear_solver = std::make_shared<LUDecomposition<Matrix, Vector>>();
                linear_solver->set_library_type("petsc");

                auto subproblem =
                    std::make_shared<utopia::MoreSorensenEigen<Matrix, Vector>>(linear_solver, eigen_solver);
                TrustRegion<Matrix, Vector> solver(subproblem);
                run_tr(this->test_functions_, solver, "TR_MS", this->verbose_);
            });
#endif  // WITH_SLEPC
#endif  // WITH_PETSC

            // this->register_experiment("PseudoTransientContinuation",
            // 	[this]() {
            // 		auto linear_solver = std::make_shared<GMRES<Matrix, Vector>>();
            // 		linear_solver->atol(1e-14);
            // 		linear_solver->max_it(10000);

            // 		PseudoContinuation<Matrix, Vector> solver(linear_solver);
            // 		solver.reset_mass_matrix(true);
            // 		run_tr(this->test_functions_, solver, "PseudoTransientContinuation", this->verbose_);
            // 	}
            // );

            // add slepcs checks
            // this->register_experiment("PseudoTR",
            // 	[this]() {
            // 		auto linear_solver = std::make_shared<GMRES<Matrix, Vector>>();
            // 		linear_solver->atol(1e-14);
            // 		linear_solver->max_it(10000);

            // 		auto eigen_solver = std::make_shared<SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL> >();
            // 		eigen_solver->solver_type("arpack");

            // 		PseudoTrustRegion<Matrix, Vector> solver(linear_solver, eigen_solver);
            // 		run_tr(this->test_functions_, solver, "PseudoTR", this->verbose_);
            // 	}
            // );

            // add slepcs checks
            // this->register_experiment("RosenbrockTR",
            // 	[this]() {
            // 		auto linear_solver = std::make_shared<GMRES<Matrix, Vector>>();
            // 		linear_solver->atol(1e-14);
            // 		linear_solver->max_it(10000);

            // 		auto eigen_solver = std::make_shared<SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL> >();
            // 		eigen_solver->solver_type("arpack");

            // 		RosenbrockTrustRegion<Matrix, Vector> solver(linear_solver, eigen_solver);
            // 		run_tr(this->test_functions_, solver, "RosenbrockTR", this->verbose_);
            // 	}
            // );

            // does not converge for most of the test cases ....
            // this->register_experiment("AffineSimilarity",
            // 	[this]() {
            // 		auto linear_solver = std::make_shared<Factorization<Matrix, Vector>>();
            // 		// linear_solver->set_type(PETSC_TAG, LU_DECOMPOSITION_TAG);  // Tags do not exist enymore s
            // 		AffineSimilarity<Matrix, Vector> solver(linear_solver);
            // 		run_tr(this->test_functions_, solver, "AffineSimilarity", this->verbose_);
            // 	}
            // );

            // 			this->register_experiment("ASTRUM",
            // 				[this]() {

            // #ifdef WITH_PETSC
            // 					auto linear_solver = std::make_shared<Factorization<Matrix, Vector>
            // >(MATSOLVERPETSC, PCLU); #else 					auto linear_solver =
            // std::make_shared<Factorization<Matrix, Vector>>(); #endif //WITH_PETSC

            // 					ASTRUM<Matrix, Vector> solver(linear_solver);
            // 					solver.reset_mass_matrix(true);
            // 					solver.verbosity_level(VERBOSITY_LEVEL_QUIET);
            // 					run_tr(this->test_functions_, solver, "ASTRUM", this->verbose_);
            // 				}
            // 			);
        }

    private:
        SizeType n_;
        std::vector<std::shared_ptr<UnconstrainedTestFunction<Matrix, Vector>>> test_functions_;
        bool verbose_{false};

        template <class TRSolver>
        static void run_tr(std::vector<std::shared_ptr<UnconstrainedTestFunction<Matrix, Vector>>> &test_functions,
                           TRSolver &solver,
                           const std::string &solv_name,
                           const bool &exp_verbose) {
            InputParameters in;
            in.set("atol", 1e-6);
            in.set("rtol", 1e-11);
            in.set("stol", 1e-14);
            in.set("stol", 1e-14);
            in.set("delta_min", 1e-13);
            in.set("max-it", 5000);
            in.set("verbose", false);

            auto params_qp = std::make_shared<InputParameters>();
            params_qp->set("atol", 1e-14);
            params_qp->set("rtol", 1e-14);
            params_qp->set("stol", 1e-14);
            auto params_qp_cast = std::static_pointer_cast<Input>(params_qp);
            in.set("linear-solver", params_qp_cast);

            solver.read(in);

            if (exp_verbose) {
                std::cout << "--------------------------------------------------------- \n";
                std::cout << "				" << solv_name << "				\n";
                std::cout << "--------------------------------------------------------- \n";
            }

            for (size_t i = 0; i < test_functions.size(); i++)
            // for(auto i =0; i < 1; i++)
            {
                Vector x_init = test_functions[i]->initial_guess();
                solver.solve(*test_functions[i], x_init);

                auto sol_status = solver.solution_status();
                // sol_status.describe(std::cout);

                const auto dim = test_functions[i]->dim();
                const auto num_its = sol_status.iterates;
                // const auto conv_reason = sol_status.reason;

                if (exp_verbose) {
                    std::cout << i << std::setw(5 - std::to_string(i).size()) << " : " << test_functions[i]->name()
                              << "_" << dim << std::right
                              << std::setw(40 - std::to_string(dim).size() - test_functions[i]->name().size())
                              << std::right << "its:  " << num_its << std::setw(5 - std::to_string(num_its).size())
                              << "  \n";

                    // if(conv_reason< 0)
                    // {
                    // 	sol_status.describe(std::cout);
                    // }
                }

                // if(test_functions[i]->exact_sol_known())
                // {
                // 	// disp(x_init);
                // 	utopia_test_assert(approxeq(x_init, test_functions[i]->exact_sol(), 1e-4));
                // }
            }
        }
    };

    static void unconstrained_opt() {
        int verbosity_level = 1;
        if (Utopia::instance().verbose()) {
            verbosity_level = 2;
        }

        if (mpi_world_size() == 1) {
#ifdef WITH_PETSC
            UnconstrainedOptimizationBenchmark<PetscMatrix, PetscVector> bench1(10);
            bench1.set_verbosity_level(verbosity_level);
            bench1.run();
#endif  // WITH_PETSC
        } else {
            std::cout << "unconstrained_opt, does not work in parallel. \n";
        }
    }

    UTOPIA_REGISTER_TEST_FUNCTION_OPTIONAL(unconstrained_opt);
}  // namespace utopia
