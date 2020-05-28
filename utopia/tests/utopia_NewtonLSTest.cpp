#include <cassert>
#include <string>
#include "utopia.hpp"
#include "utopia_Base.hpp"
#include "utopia_Benchmark.hpp"
#include "utopia_Chrono.hpp"
#include "utopia_LargeScaleIncludes.hpp"
#include "utopia_MPI.hpp"
#include "utopia_Testing.hpp"
#include "utopia_UnconstrainedBenchmark.hpp"

namespace utopia {
    template <class Matrix, class Vector>
    class NewtonLSBenchmark : public Benchmark {
    public:
        DEF_UTOPIA_SCALAR(Vector);

        std::string name() override { return "NewtonLSBenchmark benchmark."; }

        explicit NewtonLSBenchmark(const SizeType &n, const bool verbose) : n_(n), verbose_(verbose) {
            if (mpi_world_size() == 1) {
                test_functions_parallel_.resize(2);
                test_functions_parallel_[0] = std::make_shared<Poisson1D<Matrix, Vector> >(n_);
                test_functions_parallel_[1] = std::make_shared<Bratu1D<Matrix, Vector> >(n_);

                test_functions_.resize(3);
                test_functions_[0] = std::make_shared<Rosenbrock01<Matrix, Vector> >();
                test_functions_[1] = std::make_shared<QPTestFunction_2D<Matrix, Vector> >();
                test_functions_[2] = std::make_shared<Watson20<Matrix, Vector> >();
            } else {
                test_functions_parallel_.resize(2);
                test_functions_parallel_[0] = std::make_shared<Poisson1D<Matrix, Vector> >(n_ * mpi_world_size());
                test_functions_parallel_[1] = std::make_shared<Bratu1D<Matrix, Vector> >(n_ * mpi_world_size());
            }
        }

        ~NewtonLSBenchmark() override { test_functions_.clear(); }

        void initialize() override {
            // As not all backends have factorization
            // this->register_experiment("NewtonTest_FACTORIZATION",
            // 	[this]() {
            //            auto lin_solver = std::make_shared<utopia::Factorization<Matrix, Vector> >();
            //            // lin_solver->set_type(const std::string &lib, const std::string &type)

            //            Newton<Matrix, Vector> solver(lin_solver);
            //            solver.verbose(true);
            //            run_newton(this->test_functions_parallel_, solver, "NewtonTest_FACTORIZATION",
            //            this->verbose_); run_newton(this->test_functions_,  solver, "NewtonTest_FACTORIZATION",
            //            this->verbose_);
            // 	}
            // );

            this->register_experiment("NewtonTest_CG_HOMEMADE_GS", [this]() {
                auto lin_solver = std::make_shared<utopia::ConjugateGradient<Matrix, Vector, HOMEMADE> >();
                lin_solver->set_preconditioner(std::make_shared<GaussSeidel<Matrix, Vector, HOMEMADE> >());

                Newton<Matrix, Vector> solver(lin_solver);
                run_newton(this->test_functions_parallel_, solver, "NewtonTest_CG_HOMEMADE_GS", this->verbose_);
                run_newton(this->test_functions_, solver, "NewtonTest_CG_HOMEMADE_GS", this->verbose_);
            });

            this->register_experiment("NewtonTest_CG_BACKEND_GS", [this]() {
                auto lin_solver = std::make_shared<utopia::ConjugateGradient<Matrix, Vector> >();
                lin_solver->set_preconditioner(std::make_shared<GaussSeidel<Matrix, Vector> >());

                Newton<Matrix, Vector> solver(lin_solver);
                run_newton(this->test_functions_parallel_, solver, "NewtonTest_CG_BACKEND_GS", this->verbose_);
                run_newton(this->test_functions_, solver, "NewtonTest_CG_BACKEND_GS", this->verbose_);
            });

            this->register_experiment("NewtonTest_STCG_HOMEMADE_GS", [this]() {
                auto lin_solver = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();
                lin_solver->set_preconditioner(std::make_shared<GaussSeidel<Matrix, Vector, HOMEMADE> >());

                Newton<Matrix, Vector> solver(lin_solver);
                run_newton(this->test_functions_parallel_, solver, "NewtonTest_CG_HOMEMADE_GS", this->verbose_);
                run_newton(this->test_functions_, solver, "NewtonTest_CG_HOMEMADE_GS", this->verbose_);
            });

            this->register_experiment("NewtonTest_STCG_BACKEND_GS", [this]() {
                auto lin_solver = std::make_shared<utopia::SteihaugToint<Matrix, Vector> >();
                lin_solver->set_preconditioner(std::make_shared<GaussSeidel<Matrix, Vector> >());

                Newton<Matrix, Vector> solver(lin_solver);
                run_newton(this->test_functions_parallel_, solver, "NewtonTest_STCG_BACKEND_GS", this->verbose_);
                run_newton(this->test_functions_, solver, "NewtonTest_STCG_BACKEND_GS", this->verbose_);
            });

            // TODO(zulianp): : BiCGSTAB seems to have problem
            // this->register_experiment("NewtonTest_BiCGSTAB_HOMEMADE_Jacobi",
            // 	[this]() {
            //            auto lin_solver = std::make_shared<utopia::BiCGStab<Matrix, Vector, HOMEMADE> >();
            //            lin_solver->set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector> >());

            //            Newton<Matrix, Vector> solver(lin_solver);
            //            run_newton(this->test_functions_parallel_, solver, "NewtonTest_BiCGSTAB_HOMEMADE_Jacobi",
            //            this->verbose_); run_newton(this->test_functions_, solver,
            //            "NewtonTest_BiCGSTAB_HOMEMADE_Jacobi", this->verbose_);
            // 	}
            // );

            // this->register_experiment("NewtonTest_BiCGSTAB_BACKEND_Jacobi",
            // 	[this]() {
            //            auto lin_solver = std::make_shared<utopia::BiCGStab<Matrix, Vector> >();
            //            lin_solver->set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector> >());

            //            Newton<Matrix, Vector> solver(lin_solver);
            //            run_newton(this->test_functions_parallel_, solver, "NewtonTest_BiCGSTAB_BACKEND_Jacobi",
            //            this->verbose_); run_newton(this->test_functions_, solver,
            //            "NewtonTest_BiCGSTAB_BACKEND_Jacobi", this->verbose_);
            // 	}
            // );

            // working
            this->register_experiment("NewtonTest_STCG_SimpleBacktracking", [this]() {
                auto lin_solver = std::make_shared<utopia::SteihaugToint<Matrix, Vector> >();
                lin_solver->set_preconditioner(std::make_shared<GaussSeidel<Matrix, Vector> >());

                Newton<Matrix, Vector> solver(lin_solver);
                auto ls_strat = std::make_shared<utopia::SimpleBacktracking<Vector> >();
                solver.set_line_search_strategy(ls_strat);

                run_newton(
                    this->test_functions_parallel_, solver, "NewtonTest_STCG_SimpleBacktracking", this->verbose_);
                run_newton(this->test_functions_, solver, "NewtonTest_STCG_SimpleBacktracking", this->verbose_);
            });

            // working
            this->register_experiment("NewtonTest_BiCGSTAB_Backtracking", [this]() {
                auto lin_solver = std::make_shared<utopia::SteihaugToint<Matrix, Vector> >();
                lin_solver->set_preconditioner(std::make_shared<GaussSeidel<Matrix, Vector> >());

                Newton<Matrix, Vector> solver(lin_solver);
                auto ls_strat = std::make_shared<utopia::Backtracking<Vector> >();
                solver.set_line_search_strategy(ls_strat);

                run_newton(this->test_functions_parallel_, solver, "NewtonTest_BiCGSTAB_Backtracking", this->verbose_);
                run_newton(this->test_functions_, solver, "NewtonTest_BiCGSTAB_Backtracking", this->verbose_);
            });

            this->register_experiment("Inexact_Newton_test", [this]() {
                auto lin_solver = std::make_shared<utopia::SteihaugToint<Matrix, Vector> >();

                Newton<Matrix, Vector> solver(lin_solver);
                solver.forcing_strategy(InexactNewtonForcingStartegies::CAI);
                run_newton(this->test_functions_, solver, "Inexact_Newton_test", this->verbose_);
                run_newton(this->test_functions_parallel_, solver, "Inexact_Newton_test", this->verbose_);
            });
        }

    private:
        template <class Fun, class NonlinearSolver>
        static void run_newton(std::vector<std::shared_ptr<Fun> > &test_functions,
                               NonlinearSolver &solver,
                               const std::string &solv_name,
                               const bool &exp_verbose) {
            InputParameters in;
            in.set("atol", 1e-6);
            in.set("rtol", 1e-11);
            in.set("stol", 1e-14);
            in.set("stol", 1e-14);
            in.set("delta_min", 1e-13);
            in.set("max-it", 130);
            in.set("verbose", false);
            solver.read(in);

            if (exp_verbose && mpi_world_rank() == 0) {
                std::cout << "--------------------------------------------------------- \n";
                std::cout << "				" << solv_name << "				\n";
                std::cout << "--------------------------------------------------------- \n";
            }

            for (size_t i = 0; i < test_functions.size(); i++) {
                Vector x_init = test_functions[i]->initial_guess();
                solver.solve(*test_functions[i], x_init);

                auto sol_status = solver.solution_status();
                // sol_status.describe(std::cout);

                const auto dim = test_functions[i]->dim();
                const auto num_its = sol_status.iterates;
                const auto conv_reason = sol_status.reason;

                if (conv_reason <= 0) {
                    std::cerr << "[Error] NewtonLSBenchmark::run_newton(..., " << solv_name
                              << ") not converged (reason = " << ConvergenceReason::diverged_reason_string(conv_reason)
                              << " ) [[[[[[ please fix this :) ]]]]]]" << std::endl;

                    utopia_test_assert(conv_reason > 0);
                }

                // FIXME uncomment me out once fixed
                // utopia_test_assert(conv_reason > 0);

                if (exp_verbose && mpi_world_rank() == 0) {
                    std::cout << i << std::setw(5 - std::to_string(i).size()) << " : " << test_functions[i]->name()
                              << "_" << dim << std::right
                              << std::setw(60 - std::to_string(dim).size() - test_functions[i]->name().size())
                              << std::right << "its:  " << num_its << std::setw(5 - std::to_string(num_its).size())
                              << "  \n";

                    if (conv_reason < 0) {
                        sol_status.describe(std::cout);
                    }
                }

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
        std::vector<std::shared_ptr<UnconstrainedExtendedTestFunction<Matrix, Vector> > > test_functions_parallel_;
        std::vector<std::shared_ptr<UnconstrainedTestFunction<Matrix, Vector> > > test_functions_;
        SizeType n_;
        bool verbose_;
    };

    static void newton_ls() {
#ifdef WITH_PETSC
        int verbosity_level = 1;
        const int n_global = 10;
        bool alg_verbose = false;

        if (Utopia::instance().verbose()) {
            verbosity_level = 2;
        }

        NewtonLSBenchmark<PetscMatrix, PetscVector> bench1(n_global, alg_verbose);
        bench1.set_verbosity_level(verbosity_level);
        bench1.run();
#endif  // WITH_PETSC
    }

    UTOPIA_REGISTER_TEST_FUNCTION_OPTIONAL(newton_ls);
}  // namespace utopia
