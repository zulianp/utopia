#ifndef UTOPIA_SOLVER_TESTFUNCTIONND_HPP
#define UTOPIA_SOLVER_TESTFUNCTIONND_HPP

#include "utopia_TestFunctionsND.hpp"

namespace utopia {

    /**
     * @brief      Example of ND nonlinear function. Used to test nonlinear solvers.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <class Matrix, class Vector>
    class QuadraticOffsetFunction_ND final : public UnconstrainedTestFunction<Matrix, Vector> {
        using Traits = utopia::Traits<Vector>;
        using Comm = typename Traits::Communicator;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;

    public:
        // Parameter n represents the number of variables
        // Setup function f(x_0, x_1, .. x_n) = x^2 + (x + 1)^2 + ... + (x_n + n)^2
        QuadraticOffsetFunction_ND(const Comm &comm, SizeType n) : n_(n) {
            x_init_.zeros(layout(comm, Traits::decide(), n));
            x_exact_.zeros(layout(comm, Traits::decide(), n));
            /* init implemented as separate function to be able to exploit parallel region:
             * embedding it in constructor would cause compile error with cuda backend;
             * defining it as private method, would also cause compile error with cuda backend.
             */
            init_exact_solution();
        }

        void init_exact_solution() {
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
            result.identity(square_matrix_layout(layout(x_init_)), 2.0);
            return true;
        }

        Vector initial_guess() const override { return x_init_; }

        const Vector &exact_sol() const override { return x_exact_; }

        std::string name() const override { return "QuadraticOffsetFunction_ND"; }

        SizeType dim() const override { return n_; }

        bool exact_sol_known() const override { return true; }

    private:
        const SizeType n_;
        Vector x_init_;
        Vector x_exact_;
    };

}  // namespace utopia

#endif  // UTOPIA_SOLVER_TESTFUNCTIONND_HPP
