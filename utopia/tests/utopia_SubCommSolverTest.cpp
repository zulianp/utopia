#include "utopia.hpp"
#include "utopia_Newton.hpp"
#include "utopia_QPTestFunction2D.hpp"
#include "utopia_SubCommUnitTest.hpp"
#include "utopia_Testing.hpp"

namespace utopia {
    template <class Matrix, class Vector>
    class QuadraticOffsetFunction_ND final : public UnconstrainedTestFunction<Matrix, Vector> {
        using Traits = utopia::Traits<Vector>;
        using Comm = typename Traits::Communicator;
        using Scalar = typename Vector::Scalar;
        using SizeType = typename Vector::SizeType;

    public:
        // Parameter n represents the number of variables
        // Setup function f(x_0, x_1, .. x_n) = x^2 + (x - 1)^2 + ... + (x_n - n)^2
        QuadraticOffsetFunction_ND(Comm &comm, size_t n) {
            x_init_.zeros(layout(comm, Traits::decide(), n));
            x_exact_.zeros(layout(comm, Traits::decide(), n));

            const auto i_start = range(x_exact_).begin();
            auto x_view = local_view_device(x_exact_);
            parallel_for(
                local_range_device(x_exact_),
                UTOPIA_LAMBDA(const SizeType &i) { x_view.set(i, -1.0 * (i_start + i)); });
        }

        bool value(const Vector &point, Scalar &result) const override {
            const auto i_start = range(point).begin();
            auto p_view = const_local_view_device(point);

            Scalar v_sum = 0;
            parallel_reduce(
                local_range_device(point),
                UTOPIA_LAMBDA(const SizeType &i) {
                    return (p_view.get(i) + i_start + i) * (p_view.get(i) + i_start + i);
                },
                v_sum);

            result = point.comm().sum(v_sum);
            return true;
        }

        bool gradient(const Vector &point, Vector &result) const override {
            if (empty(result)) {
                result.zeros(layout(point));
            }
            const auto i_start = range(point).begin();
            auto p_view = const_local_view_device(point);
            auto r_view = local_view_device(result);
            parallel_for(
                local_range_device(point),
                UTOPIA_LAMBDA(const SizeType &i) { r_view.set(i, 2 * (p_view.get(i) + i_start + i)); });

            return true;
        }

        bool hessian(const Vector & /*point*/, Matrix &result) const override {
            result.identity(layout(x_init_.comm(), Traits::decide(), Traits::decide(), x_init_.size(), x_init_.size()),
                            2.0);
            return true;
        }

        Vector initial_guess() const override { return x_init_; }

        const Vector &exact_sol() const override { return x_exact_; }

        Scalar min_function_value() const override { return 0; }

        std::string name() const override { return "QuadraticOffsetFunction_ND"; }

        SizeType dim() const override { return 2; }

        bool exact_sol_known() const override { return true; }

    private:
        Vector x_init_;
        Vector x_exact_;
    };

    template <class Matrix, class Vector>
    class GradientDescentTest : public SubCommUnitTest<Vector> {
    public:
        void run() {
            UTOPIA_RUN_TEST(grad_descent_solve_quadratic_ND);
            UTOPIA_RUN_TEST(grad_descent_solve_quadratic_offset_ND);
        }

    private:
        using Traits = utopia::Traits<Vector>;

        void solve_and_verify(FunctionBase<Vector> &fun, Vector &x, const Vector &expected) const {
            auto solver = GradientDescent<Vector>();
            solver.dumping_parameter(0.05);
            solver.solve(fun, x);

            utopia_test_assert(approxeq(expected, x));
        }

        void grad_descent_solve_quadratic_ND() {
            constexpr size_t n = 100;
            SimpleQuadraticFunction<Matrix, Vector> fun(n);
            Vector x_init(layout(this->comm(), Traits::decide(), n));
            {
                auto x_init_view = local_view_device(x_init);
                parallel_for(
                    local_range_device(x_init), UTOPIA_LAMBDA(const SizeType &i) { x_init_view.set(i, 1.0 * i); });
            }
            Vector expected(layout(x_init), 0.0);
            solve_and_verify(fun, x_init, expected);
        }

        void grad_descent_solve_quadratic_offset_ND() {
            constexpr size_t n = 100;
            QuadraticOffsetFunction_ND<Matrix, Vector> fun(this->comm(), n);
            auto x_init = fun.initial_guess();
            solve_and_verify(fun, x_init, fun.exact_sol());
        }
    };

    template <class Matrix, class Vector>
    class NewtonTest : public SubCommUnitTest<Vector> {
    public:
        void run() {
            UTOPIA_RUN_TEST(newton_solve_quadratic_2D);
            UTOPIA_RUN_TEST(newton_solve_quadratic_offset_ND);
        }

    private:
        using Solver = utopia::LinearSolver<Matrix, Vector>;

        void solve_and_verify(UnconstrainedTestFunction<Matrix, Vector> &fun,
                              const std::shared_ptr<Solver> &linear_solver) {
            Newton<Matrix, Vector> solver(linear_solver);

            InputParameters in;
            in.set("atol", 1e-5);
            in.set("rtol", 1e-11);
            in.set("stol", 1e-14);
            in.set("stol", 1e-14);
            in.set("delta_min", 1e-13);
            in.set("max_it", 100);
            in.set("verbose", false);
            solver.read(in);

            Vector x = fun.initial_guess();
            solver.solve(fun, x);

            if (fun.exact_sol_known()) {
                utopia_test_assert(approxeq(fun.exact_sol(), x));
            } else {
                auto sol_status = solver.solution_status();
                utopia_test_assert(sol_status.reason >= 0);
            }
        }

        void newton_solve_quadratic_2D() {
            QPTestFunction_2D<Matrix, Vector> fun;

            const auto linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector>>();
            solve_and_verify(fun, linear_solver);
        }

        void newton_solve_quadratic_offset_ND() {
            constexpr size_t n = 100;
            QuadraticOffsetFunction_ND<Matrix, Vector> fun(this->comm(), n);

            const auto linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector>>();
            solve_and_verify(fun, linear_solver);
        }
    };

    void sub_comm_solver() {
        const bool verbose = Utopia::instance().verbose();
#ifdef UTOPIA_WITH_BLAS
        // Serial backend
        run_serial_test<GradientDescentTest<BlasMatrixd, BlasVectord>>();
        run_serial_test<NewtonTest<BlasMatrixd, BlasVectord>>();
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