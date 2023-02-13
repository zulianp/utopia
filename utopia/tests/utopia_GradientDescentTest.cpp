#include "utopia.hpp"
#include "utopia_AlgebraUnitTest.hpp"
#include "utopia_TestFunctionsND.hpp"
#include "utopia_Testing.hpp"

namespace utopia {
    template <class Matrix, class Vector>
    class GradientDescentTest : public AlgebraUnitTest<Vector> {
    public:
        void run() {
            UTOPIA_RUN_TEST(solve_simple_quadratic);
            UTOPIA_RUN_TEST(solve_offset_quadratic);
        }

    private:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Vector::Scalar;

        class OffsetQuadraticFunction : public Function<Matrix, Vector> {
        public:
            OffsetQuadraticFunction(const Scalar &a) : a_(a) {}

            bool value(const Vector &point, Scalar &result) const override {
                const Scalar val = norm2(point);
                // (x + a)^2 = x^2 + 2ax + a^2
                result = val * val + 2 * a_ * val + a_ * a_;
                return true;
            }

            bool gradient(const Vector &point, Vector &result) const override {
                result = 2 * (point + Vector(layout(point), a_));
                return true;
            }

            bool hessian(const Vector & /*point*/, Matrix & /*result*/) const override {
                assert(false && "Method not implemented");
                return false;
            }

        private:
            Scalar a_;
        };

        void solve_and_verify(size_t n, Function<Matrix, Vector> &fun, Scalar val_expected) {
            Vector actual(layout(this->comm(), Traits::decide(), n), 1.0);
            {
                auto a_view = local_view_device(actual);
                parallel_for(
                    local_range_device(actual), UTOPIA_LAMBDA(const SizeType &i) {
                        const auto val = a_view.get(i);
                        a_view.set(i, val * i);
                    });
            }

            auto solver = GradientDescent<Vector>();
            solver.dumping_parameter(0.05);
            solver.solve(fun, actual);

            Vector expected(layout(actual), val_expected);
            utopia_test_assert(approxeq(expected, actual));
        }

        void solve_simple_quadratic() {
            constexpr size_t n = 100;
            SimpleQuadraticFunction<Matrix, Vector> fun(n);
            solve_and_verify(n, fun, 0);
        }

        void solve_offset_quadratic() {
            constexpr size_t n = 100;
            constexpr Scalar x_offset = 2;
            OffsetQuadraticFunction fun(x_offset);
            solve_and_verify(n, fun, -x_offset);
        }
    };

    void gradient_descent() {
        const bool verbose = Utopia::instance().verbose();
#ifdef UTOPIA_WITH_BLAS
        // Serial backend
        run_serial_test<GradientDescentTest<BlasMatrixd, BlasVectord>>();
#endif  // UTOPIA_WITH_BLAS

#ifdef UTOPIA_WITH_PETSC
        run_parallel_test<GradientDescentTest<PetscMatrix, PetscVector>>(verbose);
#endif  // UTOPIA_WITH_PETSC

#ifdef UTOPIA_WITH_TRILINOS
        run_parallel_test<GradientDescentTest<TpetraMatrixd, TpetraVectord>>(verbose);
#endif  // UTOPIA_WITH_TRILINOS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(gradient_descent);

}  // namespace utopia