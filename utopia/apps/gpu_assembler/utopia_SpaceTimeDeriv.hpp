#ifndef UTOPIA_SPACE_TIME_DERIV_HPP
#define UTOPIA_SPACE_TIME_DERIV_HPP

#include "utopia_FunctionSpace.hpp"
#include "utopia_AssemblyView.hpp"

#include <string>

namespace utopia {


    template<
        class Elem,
        class Quadrature,
        class MemType = typename Elem::MemType,
        typename...
    >
    class SpaceTimeDerivView {
    public:
        static const int Dim = Elem::Dim;
        using Scalar = typename Elem::Scalar;
      //        using SizeType = typename Elem::SizeType;
        using Point      = typename Elem::Point;
        using GradValue  = typename Elem::GradValue;
        using STGradX    = typename Elem::STGradX;

        static const std::size_t NFunctions = Elem::NFunctions;
        static const std::size_t NQPoints   = Quadrature::NPoints;

        using STGradXValues        = utopia::ArrayView<STGradX, NQPoints, NFunctions>;
        using STGradPartialTValues = utopia::ArrayView<Scalar, NQPoints, NFunctions>;

        UTOPIA_INLINE_FUNCTION SpaceTimeDerivView(
            const STGradXValues &grad_x,
            const STGradXValues &grad_x_partial_t,
            const STGradPartialTValues &partial_t)
        : grad_x_(grad_x), grad_x_partial_t_(grad_x_partial_t), partial_t_(partial_t)
        {}

        UTOPIA_INLINE_FUNCTION const STGradX &grad_x(const int fun_num, const int qp_idx) const
        {
           return grad_x_(qp_idx, fun_num);
        }

        UTOPIA_INLINE_FUNCTION const STGradX &grad_x_partial_t(const int fun_num, const int qp_idx) const
        {
           return grad_x_partial_t_(qp_idx, fun_num);
        }

        UTOPIA_INLINE_FUNCTION const Scalar &partial_t(const int fun_num, const int qp_idx) const
        {
           return partial_t_(qp_idx, fun_num);
        }

        UTOPIA_INLINE_FUNCTION std::size_t n_points() const
        {
            return NQPoints;
        }

        UTOPIA_INLINE_FUNCTION constexpr static std::size_t n_functions()
        {
            return NFunctions;
        }

        UTOPIA_INLINE_FUNCTION const SpaceTimeDerivView & make(const Elem &elem) const
        {
            return *this;
        }

    private:
        const STGradXValues        &grad_x_;
        const STGradXValues        &grad_x_partial_t_;
        const STGradPartialTValues &partial_t_;
    };

    template<class FunctionSpace, class Quadrature, typename...Args>
    class SpaceTimeDeriv {
    public:
        using Elem = typename FunctionSpace::ViewDevice::Elem;

        using ViewDevice = utopia::SpaceTimeDerivView<Elem, typename Quadrature::ViewDevice>;
        using ViewHost   = utopia::SpaceTimeDerivView<Elem, typename Quadrature::ViewHost>;

        using STGradXValues        = typename ViewDevice::STGradXValues;
        using STGradPartialTValues = typename ViewDevice::STGradPartialTValues;

        static const std::size_t NFunctions = Elem::NFunctions;
        static const std::size_t NQPoints   = Quadrature::NPoints;

        SpaceTimeDeriv(const FunctionSpace &space, const Quadrature &q)
        {
            init(space, q);
        }

        ViewDevice view_device() const
        {
            return ViewDevice(grad_x_, grad_x_partial_t_, partial_t_);
        }

    private:
        STGradXValues        grad_x_;
        STGradXValues        grad_x_partial_t_;
        STGradPartialTValues partial_t_;

        void init(const FunctionSpace &space, const Quadrature &q)
        {
            Elem e;
            space.elem(0, e);

            assert(q.n_points() == NQPoints);

            for(std::size_t qp = 0; qp < NQPoints; ++qp) {
                auto p = q.point(qp);
                for(std::size_t i = 0; i < NFunctions; ++i) {
                    e.grad_x(i, q.point(qp), grad_x_(qp, i));
                    e.grad_x_partial_t(i, p, grad_x_partial_t_(qp, i));
                    partial_t_(qp, i) = e.partial_t(i, p);
                }
            }
        }
    };

    //////////////////////////////////////////////////////////////////////////////////
}


#endif //UTOPIA_SPACE_TIME_DERIV_HPP
