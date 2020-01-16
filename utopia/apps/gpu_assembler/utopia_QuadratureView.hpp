#ifndef UTOPIA_QUADRATURE_VIEW_HPP
#define UTOPIA_QUADRATURE_VIEW_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_MemType.hpp"
#include "utopia_UniformQuad4.hpp"
#include "utopia_Quadrature.hpp"

namespace utopia {

    template<
        class PointView,
        class WeightView,
        int Dim_,
        int NPoints_,
        typename...
    >
    class QuadratureView {
    public:
        using Scalar = typename Traits<WeightView>::Scalar;

        static const int Dim     = Dim_;
        static const int NPoints = NPoints_;

        UTOPIA_INLINE_FUNCTION static constexpr int n_points()
        {
            return NPoints;
        }

        UTOPIA_INLINE_FUNCTION static constexpr int dim()
        {
            return Dim;
        }

        template<class Point>
        UTOPIA_INLINE_FUNCTION void point(const int qp_idx, Point &p) const
        {
            p[0] = points_(qp_idx, 0);
            p[1] = points_(qp_idx, 1);
        }

        UTOPIA_INLINE_FUNCTION const Scalar &weight(const int qp_idx) const
        {
            return weights_[qp_idx];
        }

        QuadratureView(
            const PointView &points,
            const WeightView &weights
        ) : points_(points), weights_(weights)
        {}

    private:
        PointView  points_;
        WeightView weights_;
    };

}


#endif //UTOPIA_QUADRATURE_VIEW_HPP
