#ifndef UTOPIA_UTOPIA_FINITEDIFFERENCE_HPP
#define UTOPIA_UTOPIA_FINITEDIFFERENCE_HPP

#include <limits>
#include "utopia_Layout.hpp"
#include "utopia_Range.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    namespace internals {
        template <class Matrix, int FILL_TYPE = Traits<Matrix>::FILL_TYPE>
        class HessianFD {
        public:
            using Scalar = typename utopia::Traits<Matrix>::Scalar;

            template <class Fun, class Vector>
            static bool apply(Fun &fun, const Vector &x, const Scalar h, Matrix &H) {
                if (H.is_sparse()) {
                    return HessianFD<Matrix, FillType::SPARSE>::apply(fun, x, h, H);
                } else {
                    return HessianFD<Matrix, FillType::DENSE>::apply(fun, x, h, H);
                }
            }

            template <class Fun, class Vector>
            static bool apply_from_grad(Fun &fun, const Vector &x, const Scalar h, Matrix &H) {
                if (H.is_sparse()) {
                    return HessianFD<Matrix, FillType::SPARSE>::apply_from_grad(fun, x, h, H);
                } else {
                    return HessianFD<Matrix, FillType::DENSE>::apply_from_grad(fun, x, h, H);
                }
            }
        };

        template <class Matrix>
        class HessianFD<Matrix, FillType::DENSE> {
        public:
            using SizeType = typename Traits<Matrix>::SizeType;
            using Scalar = typename Traits<Matrix>::Scalar;

            template <class Fun, class Vector>
            static bool apply_from_grad(Fun &fun, const Vector &x, const Scalar h, Matrix &H) {
                assert(x.comm().size() == 1 && "only for serial runs");

                const SizeType n = x.size();

                auto vec_layout = layout(x);
                auto mat_layout = square_matrix_layout(vec_layout);

                H.dense(mat_layout, 0.0);
                Vector ei(vec_layout, 0.0);
                Vector g_m = ei;
                Vector g_p = ei;
                Vector x_m = ei;
                Vector x_p = ei;

                Vector h_i = ei;

                Write<Matrix> w_H(H);

                for (SizeType i = 0; i < n; ++i) {
                    {  // Scoped lock
                        const Write<Vector> ewlock(ei);
                        if (i > 0) ei.set(i - 1, 0);
                        ei.set(i, h);
                    }

                    x_m = x - ei;
                    x_p = x + ei;

                    fun.gradient(x_m, g_m);
                    fun.gradient(x_p, g_p);
                    h_i = 1. / (2 * h) * (g_p - g_m);

                    Read<Vector> r_h(h_i);
                    for (SizeType j = 0; j < n; ++j) {
                        const Scalar val = h_i.get(j);
                        H.set(i, j, val);
                    }
                }

                return true;
            }

            template <class Fun, class Vector>
            static bool apply(Fun &fun, const Vector &x, const Scalar h, Matrix &H) {
                auto vec_layout = layout(x);
                auto mat_layout = square_matrix_layout(vec_layout);

                const SizeType n = x.size();
                H.dense(mat_layout, 0.0);

                Vector ei(vec_layout, n);
                Vector ej(vec_layout, n);

                const Write<Matrix> wlock(H);

                Scalar h2 = h * h;

                const Range rr = row_range(H);
                const Range cr = col_range(H);
                const Range vr = range(x);

                for (SizeType i = 0; i < n; ++i) {
                    {  // Scoped lock
                        const Write<Vector> ewlock(ei);
                        if (i > 0 && vr.inside(i - 1)) ei.set(i - 1, 0);
                        if (vr.inside(i)) ei.set(i, h);
                    }

                    for (SizeType j = 0; j < n; ++j) {
                        {  // Scoped lock
                            const Write<Vector> ewlock(ej);

                            if (vr.inside((j - 1 + n) % n)) ej.set((j - 1 + n) % n, 0);
                            if (vr.inside(j)) ej.set(j, h);
                        }

                        Scalar fpipj, fpimj, fmipj, fmimj;
                        fun.value(x + ei + ej, fpipj);
                        fun.value(x + ei - ej, fpimj);
                        fun.value(x - ei + ej, fmipj);
                        fun.value(x - ei - ej, fmimj);

                        if (rr.inside(i) && cr.inside(j)) {
                            H.set(i, j, (fpipj - fpimj - fmipj + fmimj) / (4 * h2));
                        }
                    }
                }

                return true;
            }
        };

        template <class Matrix>
        class HessianFD<Matrix, FillType::SPARSE> {
        public:
            using Scalar = typename utopia::Traits<Matrix>::Scalar;

            template <class Fun, class Vector>
            static bool apply_from_grad(Fun &fun, const Vector &x, const Scalar h, Matrix &H) {
                assert(x.comm().size() == 1 && "only for serial runs");

                // const SizeType n = x.size();
                // Vector ei = zeros(n);
                Vector ei(layout(x), 0.0);

                Vector g_m = ei;
                Vector g_p = ei;
                Vector x_m = ei;
                Vector x_p = ei;

                Vector h_i = ei;
                h_i.read_lock();

                H *= 0.0;
                Matrix H_copy = H;

                Write<Matrix> w_H(H);

                SizeType last_i = -1;
                each_read(H_copy, [&](const SizeType &i, const SizeType &j, const Scalar &) {
                    if (last_i != i) {
                        {
                            // Scoped lock
                            Write<Vector> w(ei);
                            if (i > 0) ei.set(i - 1, 0);
                            ei.set(i, h);
                        }

                        x_m = x - ei;
                        x_p = x + ei;

                        fun.gradient(x_m, g_m);
                        fun.gradient(x_p, g_p);

                        h_i.read_unlock();
                        h_i = 1.0 / (2.0 * h) * (g_p - g_m);
                        h_i.read_lock();

                        last_i = i;
                    }

                    H.set(i, j, h_i.get(j));
                });

                h_i.read_unlock();
                return true;
            }

            template <class Fun, class Vector>
            static bool apply(Fun &fun, const Vector &x, const Scalar h, Matrix &H) {
                const SizeType n = x.size();
                assert(!empty(H) && "H has to be allocated with the sparsity pattern before calling this method");

                Vector ei(layout(x), 0.0);
                Vector ej(layout(x), 0.0);

                const Write<Matrix> wlock(H);

                const Scalar h2 = h * h;

                const Range rr = row_range(H);
                const Range cr = col_range(H);
                const Range vr = range(x);

                H *= 0.0;
                Matrix H_copy = H;

                SizeType last_i = -1;
                each_read(H_copy, [&](const SizeType &i, const SizeType &j, const Scalar &) {
                    if (last_i != i) {
                        const Write<Vector> ewlock(ei);
                        if (i > 0 && vr.inside(i - 1)) ei.set(i - 1, 0);
                        if (vr.inside(i)) ei.set(i, h);
                        last_i = i;
                    }

                    {  // Scoped lock
                        const Write<Vector> ewlock(ej);

                        if (vr.inside((j - 1 + n) % n)) ej.set((j - 1 + n) % n, 0);
                        if (vr.inside(j)) ej.set(j, h);
                    }

                    Scalar fpipj, fpimj, fmipj, fmimj;
                    fun.value(x + ei + ej, fpipj);
                    fun.value(x + ei - ej, fpimj);
                    fun.value(x - ei + ej, fmipj);
                    fun.value(x - ei - ej, fmimj);

                    H.set(i, j, (fpipj - fpimj - fmipj + fmimj) / (4 * h2));
                });

                return true;
            }
        };

    }  // namespace internals

    template <typename Scalar>
    class FiniteDifference {
    public:
        FiniteDifference(const Scalar h = 1e-5) : _h(h) {}

        template <class Fun, class Vector, class Matrix>
        bool hessian(Fun &fun, const Vector &x, Matrix &H, const bool from_grad = true) {
            if (from_grad) {
                return internals::HessianFD<Matrix>::apply_from_grad(fun, x, _h, H);
            } else {
                return internals::HessianFD<Matrix>::apply(fun, x, _h, H);
            }
        }

        template <class Fun, class Vector>
        void grad(Fun &fun, const Vector &x, Vector &g) {
            const SizeType n = x.size();
            g.values(layout(x), 0.0);
            Vector d(layout(x), 0.0);

            const Write<Vector> wlock(g);
            const Range r = range(x);

            for (SizeType i = 0; i < n; ++i) {
                {  // Scoped lock
                    Write<Vector> ewlock(d);
                    if (i > 0 && r.inside(i - 1)) d.set(i - 1, 0);
                    if (r.inside(i)) d.set(i, _h);
                }

                Scalar fpd;
                fun.value(x + d, fpd);
                Scalar fmd;
                fun.value(x - d, fmd);

                if (r.inside(i)) g.set(i, (fpd - fmd) / (2 * _h));
            }
        }

    private:
        Scalar _h;
    };
}  // namespace utopia

#endif  // UTOPIA_UTOPIA_FINITEDIFFERENCE_HPP
