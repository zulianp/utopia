#ifndef UTOPIA_SOLVER_TRIGONOMETRIC_26
#define UTOPIA_SOLVER_TRIGONOMETRIC_26

#include <cassert>
#include <vector>
#include "utopia_Function.hpp"

namespace utopia {
    template <class Matrix, class Vector>
    class Trigonometric26 final : public UnconstrainedTestFunction<Matrix, Vector> {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        Trigonometric26(const Comm &comm = Comm::get_default(), const SizeType &n_loc = 10) : n_loc_(n_loc) {
            init(comm, n_loc);
        }

        Trigonometric26(const SizeType &n_loc) : n_loc_(n_loc) { init(Comm::get_default(), n_loc); }

        void init(const Comm &comm, const SizeType &n_loc) {
            x_exact_.zeros(layout(comm, n_loc, Traits::determine()));
            auto x_layout = layout(x_exact_);

            SizeType n_global = x_layout.size(0);

            x_init_.values(x_layout, 1. / Scalar(n_global));
            x_inc_.values(x_layout, 1.0);

            {
                auto x_inc_view = view_device(x_inc_);
                parallel_for(range_device(x_inc_), UTOPIA_LAMBDA(const SizeType &i) { x_inc_view.set(i, i + 1.0); });
            }
        }

        std::string name() const override { return "Trigonometric"; }

        SizeType dim() const override { return n_loc_; }

        bool exact_sol_known() const override { return false; }

        bool value(const Vector &x, Scalar &result) const override {
            assert(local_size(x).get(0) == this->dim());

            SizeType n_global = size(x).get(0);

            Vector xcos = x;
            Vector xsin = x;

            cos_sin_x(x, xcos, xsin);

            Scalar s1 = sum(xcos);
            Vector t(layout(x), (n_global - s1));
            t = t + x_inc_ - xsin;
            t = t - e_mul(x_inc_, xcos);

            result = dot(t, t);

            return true;
        }

        static void cos_sin_x(const Vector &x, Vector &xcos, Vector &xsin) {
            auto x_view = const_view_device(x);

            auto xcos_view = view_device(xcos);
            auto xsin_view = view_device(xsin);

            parallel_for(range_device(xcos), UTOPIA_LAMBDA(const SizeType &i) {
                const auto xi = x_view.get(i);
                xcos_view.set(i, device::cos(xi));
                xsin_view.set(i, device::sin(xi));
            });
        }

        bool gradient(const Vector &x, Vector &g) const override {
            assert(local_size(x).get(0) == this->dim());

            SizeType n_global = size(x).get(0);
            Vector xcos = x;
            Vector xsin = x;

            cos_sin_x(x, xcos, xsin);

            Scalar s1 = sum(xcos);
            Vector t(layout(x), (n_global - s1));
            t = t + x_inc_ - xsin;
            t = t - e_mul(x_inc_, xcos);

            Scalar s2 = sum(t);
            g = e_mul(x_inc_, xsin) - xcos;
            g = e_mul(g, t);
            g = 2.0 * (g + (xsin * s2));

            return true;
        }

        bool hessian(const Vector &x, Matrix &H) const override {
            assert(local_size(x).get(0) == this->dim());

            if (empty(H)) {
                // H = local_values(local_size(x).get(0), local_size(x).get(0), 0.0);
                H.dense(square_matrix_layout(layout(x)), 0.0);
            }

            SizeType n_global = size(x).get(0);
            Vector xcos = x;
            Vector xsin = x;
            Vector ones(layout(x), 1.0);

            cos_sin_x(x, xcos, xsin);

            Scalar s1 = sum(xcos);
            Vector t = (n_global - s1) * ones;
            t = t + x_inc_ - xsin;
            t = t - e_mul(x_inc_, xcos);

            Scalar s2 = sum(t);

            {
                const Read<Vector> read(xsin);
                const Read<Vector> read2(xcos);
                const Read<Vector> read3(t);
                const Write<Matrix> write(H);

                H.transform_ijv(
                    [&xsin, &xcos, &t, s2, n_global](const SizeType i, const SizeType j, const Scalar &) -> Scalar {
                        Scalar val;

                        if (i == j) {
                            val = (((i + 1) * ((i + 1) + 2.)) + n_global) * xsin.get(i) * xsin.get(i);
                            val += xcos.get(i) * (xcos.get(i) - ((2. * (i + 1.)) + 2.) * xsin.get(i));
                            val += t.get(i) * (((i + 1.) * xcos.get(i)) + xsin.get(i));
                            val = 2.0 * (val + (xcos.get(i) * s2));
                        } else {
                            Scalar th = xcos.get(i);
                            Scalar xsj = xsin.get(i);
                            val = xsin.get(j) * (((n_global + (i + 1) + (j + 1)) * xsj) - th);
                            val -= (xsj * xcos.get(j));
                            val *= 2.0;
                        }

                        return val;
                    });
            }

            return true;
        }

        Vector initial_guess() const override { return x_init_; }

        const Vector &exact_sol() const override { return x_exact_; }

        Scalar min_function_value() const override { return 0; }

    private:
        SizeType n_loc_;
        Vector x_init_;
        Vector x_exact_;
        Vector x_inc_;
    };

}  // namespace utopia

#endif  // UTOPIA_SOLVER_TRIGONOMETRIC_26
