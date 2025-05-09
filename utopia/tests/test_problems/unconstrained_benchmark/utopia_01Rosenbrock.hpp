#ifndef UTOPIA_ROSENBROCK_01
#define UTOPIA_ROSENBROCK_01

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_TestFunctions.hpp"
#include "utopia_Traits.hpp"

namespace utopia {
    /**
     * @brief      Rosenbrock 2D banana function. \n
     *             The floor of the valley follows approximately the parabola \f$ y = x^2 + 1/200 \f$.
     *             The covariance matrix is not positive-definite. On the dashed line it is singular.
     *             Stepping method tend to perform at least as well as gradient methods for this function.
     *
     */
    template <class Matrix, class Vector>
    class Rosenbrock01 final : public UnconstrainedTestFunction<Matrix, Vector> {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        Rosenbrock01() {
            auto v_layout = serial_layout(dim());

            x_init_.zeros(v_layout);
            x_exact_.values(v_layout, 1.0);

            {
                const Write<Vector> write1(x_init_);
                x_init_.set(0, -1.2);
                x_init_.set(1, 1.0);
            }
        }

        bool value(const Vector &point, typename Vector::Scalar &result) const override {
            if (point.comm().size() > 1) {
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size() == this->dim());

            const Read<Vector> read(point);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);

            result = 1 + 100.0 * pow(x * x - y, 2.0) + pow(x - 1, 2.0);
            return true;
        }

        bool gradient(const Vector &point, Vector &result) const override {
            if (point.comm().size() > 1) {
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size() == this->dim());

            if (empty(result)) {
                result.zeros(layout(point));
            }

            const Read<Vector> read(point);
            const Write<Vector> write(result);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);

            result.set(0, (400.0 * x * x * x - 400 * x * y + 2.0 * x - 2.0));
            result.set(1, 200.0 * (y - x * x));
            return true;
        }

        bool hessian(const Vector &point, Matrix &result) const override {
            if (point.comm().size() > 1) {
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size() == this->dim());

            if (empty(result)) {
                result.dense(square_matrix_layout(layout(point)), 0.0);
            }

            const Read<Vector> read(point);
            const Write<Matrix> write(result);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);
            const Scalar mixed = -400.0 * x;

            result.set(0, 0, 1200 * x * x - 400 * y + 2);
            result.set(0, 1, mixed);
            result.set(1, 0, mixed);
            result.set(1, 1, 200.0);
            return true;
        }

        bool initialize_hessian(Matrix &H, Matrix &H_pre) const override {
            UTOPIA_UNUSED(H_pre);
            H.dense(square_matrix_layout(layout(x_init_)), 0.0);
            return true;
        }

        Vector initial_guess() const override { return x_init_; }

        const Vector &exact_sol() const override { return x_exact_; }

        std::string name() const override { return "Rosenbrock"; }

        SizeType dim() const override { return 2; }

        bool exact_sol_known() const override { return true; }

    private:
        Vector x_init_;
        Vector x_exact_;
    };
}  // namespace utopia
#endif  // UTOPIA_ROSENBROCK_01