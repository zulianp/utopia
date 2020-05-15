#ifndef UTOPIA_BROWN04
#define UTOPIA_BROWN04

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_Layout.hpp"
#include "utopia_TestFunctions.hpp"
#include "utopia_Traits.hpp"

namespace utopia {
    template <class Matrix, class Vector>
    class Brown04 final : public UnconstrainedTestFunction<Matrix, Vector> {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        Brown04() {
            auto v_layout = serial_layout(2);
            x_init_.values(v_layout, 1.0);
            x_exact_.zeros(v_layout);

            {
                const Write<Vector> write2(x_exact_);
                x_exact_.set(0, 1e6);
                x_exact_.set(1, 2e-6);
            }
        }

        std::string name() const override { return "Brown badly scaled"; }

        SizeType dim() const override { return 2.0; }

        bool value(const Vector &point, typename Vector::Scalar &result) const override {
            if (point.comm().size() > 1) {
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size() == 2);

            const Read<Vector> read(point);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);

            Scalar a = x - 1e6;
            Scalar b = y - 2e-6;
            Scalar c = x * y - 2.0;

            result = a * a + b * b + c * c;
            return true;
        }

        bool gradient(const Vector &point, Vector &result) const override {
            if (point.comm().size() > 1) {
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size() == 2);
            result.zeros(layout(point));

            const Read<Vector> read(point);
            const Write<Vector> write(result);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);

            Scalar a = (2.0 * x * (y * y + 1)) - (4.0 * (y + 500000));
            Scalar b = (2.0 * x * (x * y - 2.0)) + 2.0 * y - 2.0 / 500000;

            result.set(0, a);
            result.set(1, b);

            return true;
        }

        bool hessian(const Vector &point, Matrix &result) const override {
            if (point.comm().size() > 1) {
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size() == 2);
            result.dense(serial_layout(2, 2), 0.0);

            const Read<Vector> read(point);
            const Write<Matrix> write(result);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);

            const Scalar mixed = (4.0 * x * y) - 4.0;
            const Scalar a = 2.0 * (y * y + 1);
            const Scalar b = 2.0 * (x * x + 1);

            result.set(0, 0, a);
            result.set(0, 1, mixed);
            result.set(1, 0, mixed);
            result.set(1, 1, b);
            return true;
        }

        Vector initial_guess() const override { return x_init_; }

        const Vector &exact_sol() const override { return x_exact_; }

        Scalar min_function_value() const override { return 0; }

    private:
        Vector x_init_;
        Vector x_exact_;
    };
}  // namespace utopia

#endif  // UTOPIA_BROWN04