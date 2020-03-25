#ifndef UTOPIA_DUAL_ELEM_HPP
#define UTOPIA_DUAL_ELEM_HPP

#include "utopia_UniformQuad4.hpp"
#include "utopia_Edge2.hpp"

namespace utopia {

    template<class Elem>
    class DualElem {};

    template<typename Scalar, int Dim>
    class DualElem< Edge2<Scalar, Dim> > {
    public:
        using Elem = utopia::Edge2<Scalar, Dim>;
        using Point = typename Elem::Point;

        UTOPIA_INLINE_FUNCTION DualElem() {}

        UTOPIA_INLINE_FUNCTION DualElem(const Elem &elem, const Point &q)
        {
            init(elem, q);
        }

        UTOPIA_INLINE_FUNCTION void init(const Elem &elem, const Point &q)
        {
            StaticVector<Scalar, 2> buff;
            buff[0] = elem.fun(0, q);
            buff[1] = elem.fun(1, q);
            init(buff);
        }

        UTOPIA_INLINE_FUNCTION Scalar operator()(const int i) const
        {
            UTOPIA_DEVICE_ASSERT(i >= 0);
            UTOPIA_DEVICE_ASSERT(i < 2);
            return eval_[i];
        }

        template<class Evals>
        UTOPIA_INLINE_FUNCTION void init(const Evals &evals)
        {
            const Scalar f0 = evals[0];
            const Scalar f1 = evals[1];

            eval_[0] =  2 * f0 -     f1;
            eval_[1] =    - f0 + 2 * f1;

        }

    private:
        StaticVector<Scalar, 2> eval_;
    };

    template<typename Scalar>
    class DualElem< UniformQuad4<Scalar>> {
    public:
        using Elem = utopia::UniformQuad4<Scalar>;
        using Point = typename Elem::Point;

        UTOPIA_INLINE_FUNCTION DualElem() {}

        UTOPIA_INLINE_FUNCTION DualElem(const Elem &elem, const Point &q)
        {
            init(elem, q);
        }

        UTOPIA_INLINE_FUNCTION void init(const Elem &elem, const Point &q)
        {
            StaticVector<Scalar, 4> buff;
            buff[0] = elem.fun(0, q);
            buff[1] = elem.fun(1, q);
            buff[2] = elem.fun(2, q);
            buff[3] = elem.fun(3, q);
            init(buff);
        }

        UTOPIA_INLINE_FUNCTION Scalar operator()(const int i) const
        {
            UTOPIA_DEVICE_ASSERT(i >= 0);
            UTOPIA_DEVICE_ASSERT(i < 4);
            return eval_[i];
        }

        template<class Evals>
        UTOPIA_INLINE_FUNCTION void init(const Evals &evals)
        {
            const Scalar f0 = evals[0];
            const Scalar f1 = evals[1];
            const Scalar f2 = evals[2];
            const Scalar f3 = evals[3];

            eval_[0] =   4 * f0 - 2 * f1 +     f2 - 2 * f3;
            eval_[1] = - 2 * f0 + 4 * f1 - 2 * f2 +     f3;
            eval_[2] =   1 * f0 - 2 * f1 + 4 * f2 - 2 * f3;
            eval_[3] = - 2 * f0 +     f1 - 2 * f2 + 4 * f3;
        }

    private:
        StaticVector<Scalar, 4> eval_;
    };

    template<class DualElem, class Quadrature>
    class MassMatrixDual {
    public:
        static const Size_t NFunctions = DualElem::NFunctions;
        static const Size_t NQPoints   = Quadrature::NPoints;
        using Point  = typename Quadrature::Point;
        using Scalar = typename Quadrature::Scalar;

        MassMatrixDual()
        {
        }

        template<class Elem, class Dx>
        void init(const Elem &elem, DualElem &dual_elem, const Quadrature &q, const Dx &dx)
        {
            diag_.set(0.0);

            Point p;
            for(Size_t qp = 0; qp < NQPoints; ++qp)
            {
                q.point(qp, p);
                dual_elem.init(elem, p);

                for(Size_t i = 0; i < NFunctions; ++i) {
                    diag_(i) += dual_elem(i) * dual_elem(i) * dx(qp);

                    for(Size_t j = i + 1; j < NFunctions; ++j) {
                        auto v =  dual_elem(i) * dual_elem(j) * dx(qp);;
                        diag_(i) += v;
                        diag_(j) += v;
                    }
                }
            }
        }

    private:
        StaticVector<Scalar, NFunctions> diag_;
    };

}

#endif //UTOPIA_DUAL_ELEM_HPP
