#include <cassert>
#include <string>
#include "utopia.hpp"
#include "utopia_Benchmark.hpp"
#include "utopia_Chrono.hpp"
#include "utopia_ConstrainedBenchmark.hpp"
#include "utopia_MPI.hpp"
#include "utopia_Testing.hpp"
#include "utopia_assemble_laplacian_1D.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class ConstrainedOptimizationBenchmark : public Benchmark {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;

        std::string name() override { return "TR: Bound constrained optimization benchmark"; }

        // subproblem tol
        Scalar tol{1e-9};

        Scalar atol{1e-6};
        Scalar rtol{1e-8};
        Scalar stol{1e-8};

        bool enable_slow_solvers{false};

        ConstrainedOptimizationBenchmark() {
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

        ~ConstrainedOptimizationBenchmark() override { test_functions_.clear(); }

        void initialize() override {
            this->register_experiment("TR_Variable_Bound_MPRGP", [this]() {
                auto subproblem = std::make_shared<utopia::MPRGP<Matrix, Vector> >();
                subproblem->atol(tol);
                subproblem->stol(tol);
                subproblem->rtol(tol);
                subproblem->verbose(false);

                TrustRegionVariableBound<Matrix, Vector> tr_solver(subproblem);
                run_tr(this->test_functions_, tr_solver, "TR_Variable_Bound_MPRGP", this->verbose_);
            });

            this->register_experiment("TR_Variable_ProjGradient", [this]() {
                auto subproblem = std::make_shared<utopia::ProjectedGradient<Matrix, Vector> >();
                subproblem->atol(tol);
                subproblem->stol(tol);
                subproblem->rtol(tol);
                subproblem->verbose(false);

                TrustRegionVariableBound<Matrix, Vector> tr_solver(subproblem);
                run_tr(this->test_functions_, tr_solver, "TR_Variable_ProjGradient", this->verbose_);
            });

            // SemiSmooth Newton does not support dense matrices ....
            this->register_experiment("TR_Variable_SemiSmoothNewton", [this]() {
                auto lsolver = std::make_shared<LUDecomposition<Matrix, Vector> >();
                auto subproblem = std::make_shared<utopia::SemismoothNewton<Matrix, Vector> >(lsolver);
                subproblem->atol(tol);
                subproblem->stol(tol);
                subproblem->rtol(tol);
                subproblem->verbose(false);

                TrustRegionVariableBound<Matrix, Vector> tr_solver(subproblem);
                run_tr(this->test_functions_, tr_solver, "TR_Variable_SemiSmoothNewton", this->verbose_);
            });

            this->register_experiment("TR_Variable_SemiSmoothNewton", [this]() {
                auto lsolver = std::make_shared<LUDecomposition<Matrix, Vector> >();
                auto subproblem = std::make_shared<utopia::SemismoothNewton<Matrix, Vector> >(lsolver);
                subproblem->atol(tol);
                subproblem->stol(tol);
                subproblem->rtol(tol);
                subproblem->verbose(false);

                TrustRegionVariableBound<Matrix, Vector> tr_solver(subproblem);
                run_tr(this->test_functions_, tr_solver, "TR_Variable_SemiSmoothNewton", this->verbose_);
            });

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // FIXME(Alena) these are too slow for both testing and basic benchmarking
            if (enable_slow_solvers) {
                this->register_experiment("TR_Variable_PGS", [this]() {
                    auto subproblem = std::make_shared<utopia::ProjectedGaussSeidel<Matrix, Vector> >();
                    subproblem->atol(tol);
                    subproblem->stol(tol);
                    subproblem->rtol(tol);
                    subproblem->verbose(false);
                    subproblem->use_line_search(false);

                    TrustRegionVariableBound<Matrix, Vector> tr_solver(subproblem);
                    run_tr(this->test_functions_, tr_solver, "TR_Variable_PGS", this->verbose_);
                });

                this->register_experiment("TR_Variable_ProjCG", [this]() {
                    auto subproblem = std::make_shared<utopia::ProjectedConjugateGradient<Matrix, Vector> >();
                    subproblem->atol(tol);
                    subproblem->stol(tol);
                    subproblem->rtol(tol);
                    subproblem->verbose(false);

                    TrustRegionVariableBound<Matrix, Vector> tr_solver(subproblem);
                    run_tr(this->test_functions_, tr_solver, "TR_Variable_ProjCG", this->verbose_);
                });
            }

            // FIXME(Alena/Patrick) This is not converging
#ifdef UTOPIA_ENABLE_PETSC
            // this->register_experiment("TR_Variable_Tao", [this]() {
            //     auto lsolver = std::make_shared<LUDecomposition<PetscMatrix, PetscVector> >();
            //     lsolver->set_library_type("petsc");
            //     auto subproblem = std::make_shared<utopia::TaoQPSolver<PetscMatrix, PetscVector> >(lsolver);

            //     TrustRegionVariableBound<Matrix, Vector> tr_solver(subproblem);
            //     run_tr(this->test_functions_, tr_solver, "TR_Variable_Tao", this->verbose_);
            // });
#endif  // UTOPIA_ENABLE_PETSC
        }

    private:
        std::vector<std::shared_ptr<ConstrainedTestFunction<Matrix, Vector> > > test_functions_;
        bool verbose_{false};

        template <class TRSolver>
        void run_tr(std::vector<std::shared_ptr<ConstrainedTestFunction<Matrix, Vector> > > &test_functions,
                    TRSolver &solver,
                    const std::string &solv_name,
                    const bool &exp_verbose) const {
            InputParameters in;
            in.set("atol", atol);
            in.set("rtol", rtol);
            in.set("stol", stol);
            in.set("delta_min", 1e-13);
            in.set("max_it", 5000);
            in.set("verbose", false);

            // auto params_qp = std::make_shared<InputParameters>();
            // params_qp->set("atol", tol);
            // params_qp->set("rtol", tol);
            // params_qp->set("stol", tol);
            // auto params_qp_cast = std::static_pointer_cast<Input>(params_qp);
            // in.set("linear-solver", params_qp_cast);
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

                solver.set_box_constraints(test_functions[i]->box_constraints());
                solver.solve(*test_functions[i], x_init);

                auto sol_status = solver.solution_status();
                // sol_status.describe(std::cout);

                const auto dim = test_functions[i]->dim();
                const auto num_its = sol_status.iterates;
                // const auto conv_reason = sol_status.reason;

                if (exp_verbose) {
                    std::cout << i << std::setw(5 - std::to_string(i).size()) << " : " << test_functions[i]->name()
                              << "_" << dim << std::right
                              << std::setw(60 - std::to_string(dim).size() - test_functions[i]->name().size())
                              << std::right << "its:  " << num_its << std::setw(5 - std::to_string(num_its).size())
                              << "  \n";

                    // if(conv_reason< 0)
                    // {
                    // 	sol_status.describe(std::cout);
                    // }
                }

                if (test_functions[i]->exact_sol_known()) {
                    // disp(x_init);
                }
            }
        }
    };

#ifdef UTOPIA_ENABLE_PETSC
    static void constrained_opt() {
        if (mpi_world_size() == 1) {
            int verbosity_level = 1;
            if (Utopia::instance().verbose()) {
                verbosity_level = 2;
            }

            ConstrainedOptimizationBenchmark<PetscMatrix, PetscVector> bench1;
            bench1.set_verbosity_level(verbosity_level);
            bench1.run();

        } else {
            std::cout << "constrained_opt, does not work in parallel. \n";
        }
    }

    UTOPIA_REGISTER_TEST_FUNCTION_OPTIONAL(constrained_opt);

#endif  // UTOPIA_ENABLE_PETSC
}  // namespace utopia
