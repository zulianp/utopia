#ifndef UTOPIA_SOLVER_TESTFUNCTION2D_HPP
#define UTOPIA_SOLVER_TESTFUNCTION2D_HPP

#include <cassert>
#include <cmath>
#include <vector>
#include "utopia_TestFunctionsND.hpp"

namespace utopia {

    /**
     * @brief      Example of the nonlinear function. Used to test nonlinear solvers.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <class Matrix, class Vector>
    class QPTestFunction_2D final : public UnconstrainedTestFunction<Matrix, Vector> {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        QPTestFunction_2D(Comm &comm = Comm::get_default()) {
            x_init_.zeros(layout(comm, Traits::decide(), 2));
            x_exact_.zeros(layout(comm, Traits::decide(), 2));
            {
                Write<Vector> w(x_exact_);
                if (range(x_exact_).inside(0)) {
                    x_exact_.set(0, 6);
                }
                if (range(x_exact_).inside(1)) {
                    x_exact_.set(1, -7);
                }
            }
        }

        bool value(const Vector &point, Scalar &result) const override {
            const Read<Vector> read(point);
            Scalar v_sum = 0;
            if (range(point).inside(0)) {
                v_sum += (3.0 - 0.5 * point.get(0)) * (3.0 - 0.5 * point.get(0));
            }
            if (range(point).inside(1)) {
                v_sum += (point.get(1) + 7.0) * (point.get(1) + 7.0);
            }
            v_sum *= 4;
            result = point.comm().sum(v_sum);

            return true;
        }

        bool gradient(const Vector &point, Vector &result) const override {
            if (empty(result)) {
                result.zeros(layout(point));
            }
            // TODO: why this code would not work with Trilinos backend?
            if (Traits::Backend != TRILINOS) {
                const Read<Vector> read(point);
                const Write<Vector> write(result);

                if (range(point).inside(0)) {
                    result.set(0, 4.0 * (0.5 * point.get(0) - 3.0));
                }
                if (range(point).inside(1)) {
                    result.set(1, 8.0 * (point.get(1) + 7.0));
                }
            } else {
                const auto i_start = range(point).begin();
                auto p_view = const_local_view_device(point);
                auto r_view = local_view_device(result);
                parallel_for(
                    local_range_device(point), UTOPIA_LAMBDA(const SizeType &i) {
                        auto global_index = i_start + i;
                        if (global_index == 0) {
                            r_view.set(i, 4.0 * (0.5 * p_view.get(i) - 3.0));
                        } else if (global_index == 1) {
                            r_view.set(i, 8.0 * (p_view.get(i) + 7.0));
                        }
                    });
            }
            return true;
        }

        bool hessian(const Vector &point, Matrix &result) const override {
            if (empty(result)) {
                result.identity(layout(x_init_.comm(), Traits::decide(), Traits::decide(), 2, 2));
            }

            const Write<Matrix> write(result);
            if (range(point).inside(0)) {
                result.set(0, 0, 4.0);
            }
            if (range(point).inside(1)) {
                result.set(1, 1, 8.0);
            }

            return true;
        }

        Vector initial_guess() const override { return x_init_; }

        const Vector &exact_sol() const override { return x_exact_; }

        std::string name() const override { return "QPTestFunction_2D"; }

        SizeType dim() const override { return 2; }

        bool exact_sol_known() const override { return true; }

    private:
        Vector x_init_;
        Vector x_exact_;
    };

}  // namespace utopia

#endif  // UTOPIA_SOLVER_TESTFUNCTIONS2D_HPP
