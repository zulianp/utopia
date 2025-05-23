#include "utopia.hpp"
#include "utopia_SubCommUnitTest.hpp"
#include "utopia_TestProblems.hpp"
#include "utopia_Testing.hpp"

namespace utopia {

#ifdef UTOPIA_ENABLE_TRILINOS
    class TrilinosLinearSolverTest : public SubCommUnitTest<TpetraVector> {
    public:
        void run() override {
            // direct solvers (Amesos-based), no support for preconditioner yet
            UTOPIA_RUN_TEST(trilinos_factorization);
            // iterative solvers (Belos-based), test first without preconditioner
            UTOPIA_RUN_TEST(trilinos_cg);
            UTOPIA_RUN_TEST(trilinos_bicg);
            UTOPIA_RUN_TEST(trilinos_gmres);
            // then with preconditioner
#ifdef UTOPIA_ENABLE_TRILINOS_IFPACK2
            UTOPIA_RUN_TEST(trilinos_cg_ilut);
            UTOPIA_RUN_TEST(trilinos_bicg_ilut);
            UTOPIA_RUN_TEST(trilinos_gmres_ilut);
#endif
#ifdef UTOPIA_ENABLE_TRILINOS_MUELU
            UTOPIA_RUN_TEST(trilinos_cg_mg);
            UTOPIA_RUN_TEST(trilinos_bicg_mg);
            UTOPIA_RUN_TEST(trilinos_gmres_mg);
#endif  // UTOPIA_ENABLE_TRILINOS_MUELU
        }

    private:
        using Matrix = TpetraMatrix;
        using Vector = TpetraVector;

        static constexpr SizeType N = 1000;

        void trilinos_solver_test(std::shared_ptr<LinearSolver<Matrix, Vector>> solver) {
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

        void trilinos_factorization() {
            auto solver = std::make_shared<Factorization<Matrix, Vector>>();
            if (false) {
                InputParameters in;
                in.set("keep_symbolic_factorization", true);
                in.set("trans", "NOTRANS");
                in.set("verbose", false);
                solver->read(in);
            }
            trilinos_solver_test(solver);
        }

        void trilinos_cg(const std::string precond = "") {
            auto solver = std::make_shared<ConjugateGradient<Matrix, Vector>>(precond);
            {
                InputParameters in;
                in.set("block_size", 1);
                in.set("rtol", 1e-6);
                in.set("max_it", 500);
                in.set("orthogonalization", "ICGS");
                in.set("verbose", false);
                solver->read(in);
            }
            trilinos_solver_test(solver);
        }

        void trilinos_bicg(const std::string precond = "") {
            auto solver = std::make_shared<BiCGStab<Matrix, Vector>>(precond);
            {
                InputParameters in;
                in.set("rtol", 1e-6);
                in.set("max_it", 500);
                in.set("verbose", false);
                solver->read(in);
            }
            trilinos_solver_test(solver);
        }

        void trilinos_gmres(const std::string precond = "") {
            auto solver = std::make_shared<GMRES<Matrix, Vector>>(precond);
            {
                InputParameters in;
                in.set("block_size", 1);
                in.set("rtol", 1e-6);
                in.set("max_it", 500);
                in.set("max_restarts", 20);
                in.set("orthogonalization", "ICGS");
                in.set("verbose", false);
                solver->read(in);
            }
            trilinos_solver_test(solver);
        }

        void trilinos_minres(const std::string precond = "") {
            auto solver = std::make_shared<MINRES<Matrix, Vector>>(precond);
            {
                InputParameters in;
                in.set("block_size", 1);
                in.set("rtol", 1e-6);
                in.set("max_it", 500);
                in.set("verbose", false);
                solver->read(in);
            }
            trilinos_solver_test(solver);
        }

#ifdef UTOPIA_ENABLE_TRILINOS_IFPACK2
        void trilinos_cg_ilut() { trilinos_cg("ILUT"); }
        void trilinos_bicg_ilut() { trilinos_bicg("ILUT"); }
        void trilinos_gmres_ilut() { trilinos_gmres("ILUT"); }
        void trilinos_minres_ilut() { trilinos_minres("ILUT"); }
#endif  // UTOPIA_ENABLE_TRILINOS_IFPACK2

#ifdef UTOPIA_ENABLE_TRILINOS_MUELU
        void trilinos_cg_mg() { trilinos_cg("MueLu"); }
        void trilinos_bicg_mg() { trilinos_bicg("MueLu"); }
        void trilinos_gmres_mg() { trilinos_gmres("MueLu"); }
        void trilinos_minres_mg() { trilinos_minres("MueLu"); }
#endif  // UTOPIA_ENABLE_TRILINOS_MUELU
    };

    void trilinos_solver() {
        const bool verbose = Utopia::instance().verbose();
        run_parallel_test<TrilinosLinearSolverTest>(verbose);
    }

    UTOPIA_REGISTER_TEST_FUNCTION(trilinos_solver);

#endif  // UTOPIA_ENABLE_TRILINOS
}  // namespace utopia