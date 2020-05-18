#ifndef UTOPIA_QUADRATURE_VIEW_HPP
#define UTOPIA_QUADRATURE_VIEW_HPP

#include "utopia_Accessor.hpp"
#include "utopia_Base.hpp"
#include "utopia_MemType.hpp"
#include "utopia_Quadrature.hpp"
#include "utopia_Traits.hpp"
#include "utopia_UniformQuad4.hpp"

namespace utopia {

    template <class PointView, class WeightView, int Dim_, int NPoints_, typename...>
    class QuadratureView {
    public:
        using Scalar = typename Traits<WeightView>::Scalar;
        using Point = utopia::StaticVector<Scalar, Dim_>;

        static const int Dim = Dim_;
        static const int NPoints = NPoints_;

        using PA = utopia::Accessor<PointView>;

        UTOPIA_INLINE_FUNCTION static constexpr int n_points() { return NPoints; }

        UTOPIA_INLINE_FUNCTION static constexpr int dim() { return Dim; }

        // template<class Point>
        // UTOPIA_INLINE_FUNCTION void point(const int qp_idx, Point &p) const
        // {
        //     for(int d = 0; d < Dim; ++d) {
        //         p[d] = PA::get(points_, qp_idx, d);
        //     }
        // }

        UTOPIA_INLINE_FUNCTION const Point &point(const int qp_idx) const {
            // Point p;
            // point(qp_idx, p);
            // return p;
            return points_[qp_idx];
        }

        UTOPIA_INLINE_FUNCTION const Scalar &weight(const int qp_idx) const { return weights_[qp_idx]; }

        UTOPIA_INLINE_FUNCTION QuadratureView() = default;

        UTOPIA_INLINE_FUNCTION QuadratureView(const PointView &points, const WeightView &weights)
            : points_(points), weights_(weights) {}

        UTOPIA_INLINE_FUNCTION PointView &points() { return points_; }
        UTOPIA_INLINE_FUNCTION WeightView &weights() { return weights_; }

        UTOPIA_INLINE_FUNCTION const PointView &points() const { return points_; }
        UTOPIA_INLINE_FUNCTION const WeightView &weights() const { return weights_; }

    private:
        PointView points_;
        WeightView weights_;
    };

}  // namespace utopia

#endif  // UTOPIA_QUADRATURE_VIEW_HPP
