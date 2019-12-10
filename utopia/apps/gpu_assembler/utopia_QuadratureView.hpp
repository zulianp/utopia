#ifndef UTOPIA_QUADRATURE_VIEW_HPP
#define UTOPIA_QUADRATURE_VIEW_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_kokkos_Traits.hpp"
#include "utopia_MemType.hpp"
#include "utopia_UniformQuad4.hpp"

namespace utopia {

    template<typename T, int Order, int Dim = T::Dim, typename ...>
    class Quadrature {};

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

    template<typename Scalar_>
    class Quadrature<UniformQuad4<Scalar_>, 2, 2> {
    public:
        using Scalar = Scalar_;
        using ExecutionSpace = TpetraVector::vector_type::execution_space;

        //Host-mirror types
        using DualPointView   = Kokkos::DualView<Scalar **, ExecutionSpace>;
        using DualWeightView  = Kokkos::DualView<Scalar *, ExecutionSpace>;

        //Device types TODO
        using PointView   = Kokkos::View<Scalar **, ExecutionSpace>;
        using WeightView  = Kokkos::View<Scalar *, ExecutionSpace>;


        static const int Order   = 2;
        static const int Dim     = 2;
        static const int NPoints = 6;

        using DeviceView  = utopia::QuadratureView<PointView, WeightView, Dim, NPoints>;

        void init()
        {
            auto host_points  = points_.view_host();
            auto host_weights = weights_.view_host();

            host_points(0, 0) = 0.5;
            host_points(0, 1) = 0.5;
            host_points(1, 0) = 0.9304589153964795245728880523899,
            host_points(1, 1) = 0.5;
            host_points(2, 0) = 0.72780186391809642112479237299488;
            host_points(2, 1) = 0.074042673347699754349082179816666;
            host_points(3, 0) = 0.72780186391809642112479237299488;
            host_points(3, 1) = 0.92595732665230024565091782018333;
            host_points(4, 0) = 0.13418502421343273531598225407969;
            host_points(4, 1) = 0.18454360551162298687829339850317;
            host_points(5, 0) = 0.13418502421343273531598225407969;
            host_points(5, 1) = 0.81545639448837701312170660149683;

            host_weights(0) = 0.28571428571428571428571428571428;
            host_weights(1) = 0.10989010989010989010989010989011;
            host_weights(2) = 0.14151805175188302631601261486295;
            host_weights(3) = 0.14151805175188302631601261486295;
            host_weights(4) = 0.16067975044591917148618518733485;
            host_weights(5) = 0.16067975044591917148618518733485;
        }

        Quadrature()
        : points_("Quadrature::points", NPoints, Dim),
          weights_("Quadrature::weights", NPoints)
        {
            init();
        }

        inline DeviceView view_device() const
        {
            return DeviceView(points_.view_device(), weights_.view_device());
        }

    private:
        DualPointView points_;
        DualWeightView weights_;
    };
}


#endif //UTOPIA_QUADRATURE_VIEW_HPP
