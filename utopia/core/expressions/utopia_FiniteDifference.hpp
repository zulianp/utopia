//
// Created by Patrick Zulian on 29/05/15.
//

#ifndef UTOPIA_UTOPIA_FINITEDIFFERENCE_HPP
#define UTOPIA_UTOPIA_FINITEDIFFERENCE_HPP

#include <limits>
#include "utopia_Traits.hpp"

namespace utopia {

    namespace internals {
        template<class Matrix, int FILL_TYPE = Matrix::SPARSE>
        class HessianFD { };


        template<class Matrix>
        class HessianFD<Matrix, FillType::DENSE> {
        public:
            DEF_UTOPIA_SCALAR(Matrix);

            template<class Fun, class Vector>
            void apply(Fun &fun, const Vector &x, const Scalar h, Matrix &H) {
                auto n = x.size().get(0);
                H = zeros(n, n);
                Vector ei = zeros(n);
                Vector ej = zeros(n);

                const Write <Matrix> wlock(H);

                Scalar h2 = h * h;

                const Range rr = rowRange(H);
                const Range cr = colRange(H);
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
            }
        };


        template<class Matrix>
        class HessianFD<Matrix, FillType::SPARSE> {
        public:
            DEF_UTOPIA_SCALAR(Matrix);

            template<class Fun, class Vector>
            void apply(Fun & /*fun*/, const Vector & /*x*/, const Scalar  /*h*/, Matrix & /*H */) {
                assert(false); //TODO implement me
//
//                if(is_empty(H)) {
//                    //....
//                } else {
//
//                }
            }
        };
    }

    template<typename Scalar>
    class FiniteDifference {
    public:

        FiniteDifference(const Scalar h = 1e-5)
                : _h(h) { }

        template<class Fun, class Vector, class Matrix>
        void hessian(Fun &fun, const Vector &x, Matrix &H) {
            internals::HessianFD<Matrix> diff;
            diff.apply(fun, x, _h, H);
        }


        template<class Fun, class Vector>
        void grad(Fun &fun, const Vector &x, Vector &g) {
            auto n = x.size().get(0);
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
