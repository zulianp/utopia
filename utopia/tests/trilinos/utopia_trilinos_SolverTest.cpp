#include "utopia.hpp"
#include "utopia_SubCommUnitTest.hpp"
#include "utopia_TestProblems.hpp"
#include "utopia_Testing.hpp"

namespace utopia {

#ifdef UTOPIA_WITH_TRILINOS
    class TrilinosLinearSolverTest : public SubCommUnitTest<TpetraVector> {
    public:
        void run() {
            UTOPIA_RUN_TEST(trilinos_cg);
            UTOPIA_RUN_TEST(trilinos_gmres);
#ifdef UTOPIA_WITH_TRILINOS_IFPACK2
            UTOPIA_RUN_TEST(trilinos_cg_ilut);
            UTOPIA_RUN_TEST(trilinos_gmres_ilut);
#endif
#ifdef UTOPIA_WITH_TRILINOS_MUELU
            UTOPIA_RUN_TEST(trilinos_cg_mg);
            UTOPIA_RUN_TEST(trilinos_gmres_mg);
#endif  // UTOPIA_WITH_TRILINOS_MUELU
        }

    private:
        using Matrix = TpetraMatrix;
        using Vector = TpetraVector;
        using Traits = Traits<Vector>;

        static constexpr SizeType N = 1000;

        void belos_solver_test(std::shared_ptr<BelosSolver<Matrix, Vector>> solver, const std::string &precond) {
            if (precond.empty()) {
                utopia_test_assert(solver->get_preconditioner_name().compare(UTOPIA_PCNONE) == 0);
            } else {
                utopia_test_assert(solver->get_preconditioner_name().compare(precond) == 0);
            }

            Poisson1D<Matrix, Vector> fun(N, 2);
            Vector x = fun.initial_guess();
            Vector rhs;
            fun.get_rhs(rhs);
            Matrix A;
            fun.hessian(x, A);
            rhs *= 0.00001;
            x = rhs;

            solver->update(std::make_shared<Matrix>(A));
            solver->apply(rhs, x);

            utopia_test_assert(approxeq(rhs, A * x, 1e-5));
        }

        void trilinos_cg(const std::string precond = "") {
            {
                // 1st test case: specify preconditioner in solver constructor
                auto solver = std::make_shared<ConjugateGradient<Matrix, Vector>>(precond);

                InputParameters in;
                in.set("block_size", 1);
                in.set("rtol", 1e-6);
                in.set("max_it", 500);
                in.set("orthogonalization", "ICGS");
                in.set("verbose", false);
                solver->read(in);

                belos_solver_test(solver, precond);
            }
            {
                // 2nd test case: specify preconditioner as config parameter
                auto solver = std::make_shared<ConjugateGradient<Matrix, Vector>>();

                InputParameters in;
                in.set("block_size", 1);
                in.set("rtol", 1e-6);
                in.set("max_it", 500);
                in.set("orthogonalization", "ICGS");
                in.set("verbose", false);
                in.set("pc_type", precond);
                solver->read(in);

                belos_solver_test(solver, precond);
            }
        }

        void trilinos_gmres(const std::string precond = "") {
            {
                // 1st test case: specify preconditioner in solver constructor
                auto solver = std::make_shared<GMRES<Matrix, Vector>>(precond);

                InputParameters in;
                in.set("block_size", 1);
                in.set("rtol", 1e-6);
                in.set("max_it", 500);
                in.set("max_restarts", 20);
                in.set("orthogonalization", "ICGS");
                in.set("verbose", false);
                solver->read(in);

                belos_solver_test(solver, precond);
            }
            {
                // 2nd test case: specify preconditioner as config parameter
                auto solver = std::make_shared<GMRES<Matrix, Vector>>();

                InputParameters in;
                in.set("block_size", 1);
                in.set("rtol", 1e-6);
                in.set("max_it", 500);
                in.set("max_restarts", 20);
                in.set("orthogonalization", "ICGS");
                in.set("verbose", false);
                in.set("pc_type", precond);
                solver->read(in);

                belos_solver_test(solver, precond);
            }
        }

        void trilinos_minres(const std::string precond = "") {
            {
                // 1st test case: specify preconditioner in solver constructor
                auto solver = std::make_shared<MINRES<Matrix, Vector>>(precond);

                InputParameters in;
                in.set("block_size", 1);
                in.set("rtol", 1e-6);
                in.set("max_it", 500);
                in.set("verbose", false);
                solver->read(in);

                belos_solver_test(solver, precond);
            }
            {
                // 2nd test case: specify preconditioner as config parameter
                auto solver = std::make_shared<MINRES<Matrix, Vector>>();

                InputParameters in;
                in.set("block_size", 1);
                in.set("rtol", 1e-6);
                in.set("max_it", 500);
                in.set("verbose", false);
                in.set("pc_type", precond);
                solver->read(in);

                belos_solver_test(solver, precond);
            }
        }

#ifdef UTOPIA_WITH_TRILINOS_IFPACK2
        void trilinos_cg_ilut() { trilinos_cg(UTOPIA_PCILUT); }
        void trilinos_gmres_ilut() { trilinos_gmres(UTOPIA_PCILUT); }
        void trilinos_minres_ilut() { trilinos_minres(UTOPIA_PCILUT); }
#endif  // UTOPIA_WITH_TRILINOS_IFPACK2

#ifdef UTOPIA_WITH_TRILINOS_MUELU
        void trilinos_cg_mg() { trilinos_cg(UTOPIA_PCMULTIGRID); }
        void trilinos_gmres_mg() { trilinos_gmres(UTOPIA_PCMULTIGRID); }
        void trilinos_minres_mg() { trilinos_minres(UTOPIA_PCMULTIGRID); }
#endif  // UTOPIA_WITH_TRILINOS_MUELU
    };

    void trilinos_solver() {
        const bool verbose = Utopia::instance().verbose();
        run_parallel_test<TrilinosLinearSolverTest>(verbose);
    }

    UTOPIA_REGISTER_TEST_FUNCTION(trilinos_solver);

#endif  // UTOPIA_WITH_TRILINOS
}  // namespace utopia