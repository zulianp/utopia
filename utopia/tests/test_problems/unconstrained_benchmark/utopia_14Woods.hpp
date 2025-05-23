#ifndef UTOPIA_SOLVER_WOODS_TESTFUNCTION_HPP
#define UTOPIA_SOLVER_WOODS_TESTFUNCTION_HPP

#include <cassert>
#include <vector>
#include "utopia_Function.hpp"

namespace utopia {
    template <class Matrix, class Vector>
    class Woods14 final : public UnconstrainedTestFunction<Matrix, Vector> {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        Woods14() {
            auto v_layout = serial_layout(dim());

            x_init_.zeros(v_layout);
            x_exact_.zeros(v_layout);

            {
                const Write<Vector> write1(x_init_);
                const Write<Vector> write2(x_exact_);

                x_init_.set(0, -3.0);
                x_init_.set(1, -1.0);
                x_init_.set(2, -3.0);
                x_init_.set(3, -1.0);

                x_exact_.set(0, 1.0);
                x_exact_.set(1, 1.0);
                x_exact_.set(2, 1.0);
                x_exact_.set(3, 1.0);
            }
        }

        std::string name() const override { return "Wood"; }

        SizeType dim() const override { return 4; }

        bool value(const Vector &point, typename Vector::Scalar &result) const override {
            if (point.comm().size() > 1) {
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size() == 4);

            {
                const Read<Vector> read(point);

                const Scalar w = point.get(0);
                const Scalar x = point.get(1);
                const Scalar y = point.get(2);
                const Scalar z = point.get(3);

                result = 100 * pow(x - w * w, 2) + pow(w - 1, 2) + 90 * pow(z - y * y, 2) + pow(1 - y, 2) +
                         10.1 * (pow(x - 1, 2) + pow(z - 1, 2)) + 19.8 * (x - 1) * (z - 1);
            }
            return true;
        }

        bool gradient(const Vector &point, Vector &result) const override {
            if (point.comm().size() > 1) {
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size() == 4);
            result.zeros(layout(point));

            {
                const Read<Vector> read(point);
                const Write<Vector> write(result);

                const Scalar w = point.get(0);
                const Scalar x = point.get(1);
                const Scalar y = point.get(2);
                const Scalar z = point.get(3);

                result.set(0, 2 * w - 400 * w * (-w * w + x) - 2);
                result.set(1, -200 * w * w + (1101 * x) / 5 + (99 * z) / 5 - 40);
                result.set(2, 2 * y - 360 * y * (-y * y + z) - 2);
                result.set(3, -180 * y * y + (99 * x) / 5 + (1001 * z) / 5 - 40);
            }
            return true;
        }

        bool hessian(const Vector &point, Matrix &result) const override {
            if (point.comm().size() > 1) {
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size() == 4);
            result.dense(serial_layout(4, 4));

            const Read<Vector> read(point);
            const Write<Matrix> write(result);

            const Scalar w = point.get(0);
            const Scalar x = point.get(1);
            const Scalar y = point.get(2);
            const Scalar z = point.get(3);

            result.set(0, 0, 1200 * w * w - 400 * x + 2);
            result.set(0, 1, -400 * w);
            result.set(0, 2, 0);
            result.set(0, 3, 0);

            result.set(1, 0, -400 * w);
            result.set(1, 1, 1101 / 5);
            result.set(1, 2, 0);
            result.set(1, 3, 99 / 5);

            result.set(2, 0, 0);
            result.set(2, 1, 0);
            result.set(2, 2, 1080 * y * y - 360 * z + 2);
            result.set(2, 3, -360 * y);

            result.set(3, 0, 0);
            result.set(3, 1, 99 / 5);
            result.set(3, 2, -360 * y);
            result.set(3, 3, 1001 / 5);

            return true;
        }

        Vector initial_guess() const override { return x_init_; }

        const Vector &exact_sol() const override { return x_exact_; }

    private:
        Vector x_init_;
        Vector x_exact_;
    };

}  // namespace utopia

#endif  // UTOPIA_SOLVER_WOODS_TESTFUNCTION_HPP
