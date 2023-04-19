#include "utopia.hpp"
#include "utopia_Newton.hpp"
#include "utopia_QPTestFunction2D.hpp"
#include "utopia_QPTestFunctionND.hpp"
#include "utopia_SubCommUnitTest.hpp"
#include "utopia_Testing.hpp"

namespace utopia {
    template <class Matrix, class Vector>
    class GradientDescentTest : public SubCommUnitTest<Vector> {
    public:
        void run() override {
            UTOPIA_RUN_TEST(grad_descent_solve_quadratic_2D);
            UTOPIA_RUN_TEST(grad_descent_solve_quadratic_ND);
        }

    private:
        void solve_and_verify(UnconstrainedTestFunction<Matrix, Vector> &fun) const {
            GradientDescent<Vector> solver;
            solver.damping_parameter(0.05);

            InputParameters in;
            in.set("atol", 1e-8);
            solver.read(in);

            Vector x = fun.initial_guess();
            solver.solve(fun, x);

            utopia_assert(fun.exact_sol_known());
            utopia_test_assert(approxeq(fun.exact_sol(), x));
        }

        void grad_descent_solve_quadratic_2D() {
            QPTestFunction_2D<Matrix, Vector> fun(this->comm());
            solve_and_verify(fun);
        }

        void grad_descent_solve_quadratic_ND() {
            constexpr SizeType n = 100;
            QuadraticOffsetFunction_ND<Matrix, Vector> fun(this->comm(), n);
            solve_and_verify(fun);
        }
    };

    template <class Matrix, class Vector>
    class NewtonTest : public SubCommUnitTest<Vector> {
    public:
        void run() {
            UTOPIA_RUN_TEST(newton_solve_quadratic_2D);
            UTOPIA_RUN_TEST(newton_solve_quadratic_ND);
        }

    private:
        using Solver = utopia::LinearSolver<Matrix, Vector>;

        static constexpr SizeType N = 100;

        void solve_and_verify(UnconstrainedTestFunction<Matrix, Vector> &fun,
                              const std::shared_ptr<Solver> &linear_solver) const {
            Newton<Matrix, Vector> solver(linear_solver);

            InputParameters in;
            in.set("atol", 1e-8);
            in.set("verbose", false);
            solver.read(in);

            Vector x = fun.initial_guess();
            solver.solve(fun, x);

            utopia_assert(fun.exact_sol_known());
            utopia_test_assert(approxeq(fun.exact_sol(), x));
        }

        void newton_solve_quadratic_2D() const {
            QPTestFunction_2D<Matrix, Vector> fun(this->comm());
            solve_and_verify(fun, std::make_shared<ConjugateGradient<Matrix, Vector>>());
            solve_and_verify(fun, std::make_shared<GMRES<Matrix, Vector>>());
        }

        void newton_solve_quadratic_ND() const {
            QuadraticOffsetFunction_ND<Matrix, Vector> fun(this->comm(), N);
            solve_and_verify(fun, std::make_shared<ConjugateGradient<Matrix, Vector>>());
            solve_and_verify(fun, std::make_shared<GMRES<Matrix, Vector>>());
        }
    };

    void sub_comm_solver() {
        const bool verbose = Utopia::instance().verbose();
#ifdef UTOPIA_WITH_BLAS
        // Serial backend
        run_serial_test<GradientDescentTest<BlasMatrixd, BlasVectord>>();
        // Disable NewtonTest beacuse Blas backend does not support GMRES linear solver.
        // run_serial_test<NewtonTest<BlasMatrixd, BlasVectord>>();
#endif  // UTOPIA_WITH_BLAS

#ifdef UTOPIA_WITH_PETSC
        run_parallel_test<GradientDescentTest<PetscMatrix, PetscVector>>(verbose);
        run_parallel_test<NewtonTest<PetscMatrix, PetscVector>>(verbose);
#endif  // UTOPIA_WITH_PETSC

#ifdef UTOPIA_WITH_TRILINOS
        run_parallel_test<GradientDescentTest<TpetraMatrixd, TpetraVectord>>(verbose);
        run_parallel_test<NewtonTest<TpetraMatrix, TpetraVector>>(verbose);
#endif  // UTOPIA_WITH_TRILINOS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(sub_comm_solver);

}  // namespace utopia