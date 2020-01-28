#ifndef UTOPIA_UTOPIA_FINITEDIFFERENCE_HPP
#define UTOPIA_UTOPIA_FINITEDIFFERENCE_HPP

#include <limits>
#include "utopia_Traits.hpp"

namespace utopia {

    namespace internals {
        template<class Matrix, int FILL_TYPE = Traits<Matrix>::FILL_TYPE>
        class HessianFD {
        public:
            DEF_UTOPIA_SCALAR(Matrix);

            template<class Fun, class Vector>
            static bool apply(Fun &fun, const Vector &x, const Scalar h, Matrix &H) {

                if(H.is_sparse()) {
                    return HessianFD<Matrix, FillType::SPARSE>::apply(fun, x, h, H);
                } else {
                    return HessianFD<Matrix, FillType::DENSE>::apply(fun, x, h, H);
                }
            }

        };


        template<class Matrix>
        class HessianFD<Matrix, FillType::DENSE> {
        public:
            DEF_UTOPIA_SCALAR(Matrix);

            template<class Fun, class Vector>
            static bool apply(Fun &fun, const Vector &x, const Scalar h, Matrix &H) {
                auto n = x.size();
                H = zeros(n, n);
                Vector ei = zeros(n);
                Vector ej = zeros(n);

                const Write <Matrix> wlock(H);

                Scalar h2 = h * h;

                const Range rr = row_range(H);
                const Range cr = col_range(H);
                const Range vr = range(x);

                for (auto i = 0; i < n; ++i) {
                    { //Scoped lock
                        const Write <Vector> ewlock(ei);
                        if (i > 0 && vr.inside(i-1)) ei.set(i - 1, 0);
                        if(vr.inside(i))             ei.set(i, h);
                    }

                    for (auto j = 0; j < n; ++j) {
                        { //Scoped lock
                            const Write <Vector> ewlock(ej);

                            if(vr.inside((j - 1 + n) % n)) ej.set((j - 1 + n) % n, 0);
                            if(vr.inside(j))               ej.set(j, h);
                        }

                        Scalar fpipj, fpimj, fmipj, fmimj;
                        fun.value(x + ei + ej, fpipj);
                        fun.value(x + ei - ej, fpimj);
                        fun.value(x - ei + ej, fmipj);
                        fun.value(x - ei - ej, fmimj);

                        if(rr.inside(i) && cr.inside(j)) {
                            H.set(i, j, (   fpipj
                                            - fpimj
                                            - fmipj
                                            + fmimj
                                        ) / (4 * h2));
                        }
                    }
                }

                return true;
            }
        };


        template<class Matrix>
        class HessianFD<Matrix, FillType::SPARSE> {
        public:
            DEF_UTOPIA_SCALAR(Matrix);

            template<class Fun, class Vector>
            static bool apply(Fun &fun, const Vector &x, const Scalar h, Matrix &H) {
                auto n = x.size();
                assert(!empty(H) && "H has to be allocated with the sparsity pattern before calling this method");
                Vector ei = zeros(n);
                Vector ej = zeros(n);

                const Write <Matrix> wlock(H);

                Scalar h2 = h * h;

                const Range rr = row_range(H);
                const Range cr = col_range(H);
                const Range vr = range(x);

                H *= 0.0;
                Matrix H_copy = H;

                SizeType last_i = -1;
                each_read(H_copy, [&](const SizeType &i, const SizeType &j, const Scalar &) {
                    if(last_i != i) {
                        const Write <Vector> ewlock(ei);
                        if (i > 0 && vr.inside(i-1)) ei.set(i - 1, 0);
                        if(vr.inside(i))             ei.set(i, h);
                        last_i = i;
                    }

                    { //Scoped lock
                        const Write <Vector> ewlock(ej);

                        if(vr.inside((j - 1 + n) % n)) ej.set((j - 1 + n) % n, 0);
                        if(vr.inside(j))               ej.set(j, h);
                    }

                    Scalar fpipj, fpimj, fmipj, fmimj;
                    fun.value(x + ei + ej, fpipj);
                    fun.value(x + ei - ej, fpimj);
                    fun.value(x - ei + ej, fmipj);
                    fun.value(x - ei - ej, fmimj);

                    H.set(i, j, (fpipj
                        - fpimj
                        - fmipj
                        + fmimj
                        ) / (4 * h2));
                });

                return true;
            }
        };


        // template<class Matrix>
        // class HessianFD<Matrix, FillType::SPARSE> {
        // public:
        //     DEF_UTOPIA_SCALAR(Matrix);

        //     template<class Fun, class Vector>
        //     static bool apply(Fun & /*fun*/, const Vector & /*x*/, const Scalar  /*h*/, Matrix & /*H */) {
        //         std::cerr << ("[Warning] HessianFD not implemented for sparse matrices") << std::endl;
        //         return false;

        //     }
        // };
    }

    template<typename Scalar>
    class FiniteDifference {
    public:

        FiniteDifference(const Scalar h = 1e-5)
                : _h(h) { }

        template<class Fun, class Vector, class Matrix>
        bool hessian(Fun &fun, const Vector &x, Matrix &H) {
            return internals::HessianFD<Matrix>::apply(fun, x, _h, H);
        }


        template<class Fun, class Vector>
        void grad(Fun &fun, const Vector &x, Vector &g) {
            auto n = x.size();
            g = values(n, 1, 0.0);
            Vector d = zeros(n);

            const Write <Vector> wlock(g);
            const Range r = range(x);

            for (auto i = 0; i < n; ++i) {
                { //Scoped lock
                    Write <Vector> ewlock(d);
                    if (i > 0 && r.inside(i-1)) d.set(i - 1, 0);
                    if(r.inside(i))             d.set(i, _h);
                }

                Scalar fpd;
                fun.value(x + d, fpd);
                Scalar fmd;
                fun.value(x - d, fmd);

                if(r.inside(i)) g.set(i, (fpd - fmd) / (2 * _h));
            }
        }

    private:
        Scalar _h;

    };
}

#endif //UTOPIA_UTOPIA_FINITEDIFFERENCE_HPP
