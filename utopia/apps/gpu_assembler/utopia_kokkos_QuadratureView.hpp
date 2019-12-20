#ifndef UTOPIA_KOKKOS_QUADRATURE_VIEW_HPP
#define UTOPIA_KOKKOS_QUADRATURE_VIEW_HPP

#include "utopia_QuadratureView.hpp"
#include "utopia_kokkos_Traits.hpp"

namespace utopia {

    template<typename Scalar_>
    class Quadrature<UniformQuad4<Scalar_>, 2, 2> {
    public:
        using Scalar = Scalar_;
        using ExecutionSpace = TpetraVector::vector_type::execution_space;


        //Host-mirror types
        using DualPointView   = Kokkos::DualView<Scalar **, ExecutionSpace>;
        using DualWeightView  = Kokkos::DualView<Scalar *, ExecutionSpace>;

        using HostSpace = typename DualPointView::t_host;

        //Device types TODO
        using PointView   = Kokkos::View<Scalar **, ExecutionSpace>;
        using WeightView  = Kokkos::View<Scalar *, ExecutionSpace>;

        using PointViewHost   = Kokkos::View<Scalar **, HostSpace>;
        using WeightViewHost  = Kokkos::View<Scalar *, HostSpace>;


        static const int Order   = 2;
        static const int Dim     = 2;
        static const int NPoints = 6;

        using ViewDevice  = utopia::QuadratureView<PointView, WeightView, Dim, NPoints>;
        using ViewHost    = utopia::QuadratureView<PointViewHost, WeightViewHost, Dim, NPoints>;

        void init()
        {
            auto host_points  = points_.view_host();
            auto host_weights = weights_.view_host();
            Quad4Quadrature<Scalar, Order, Dim, NPoints>::get(host_points, host_weights);
        }

        Quadrature()
        : points_("Quadrature::points", NPoints, Dim),
          weights_("Quadrature::weights", NPoints)
        {
            init();
        }

        inline ViewDevice view_device() const
        {
            return ViewDevice(points_.view_device(), weights_.view_device());
        }

        inline ViewHost view_host() const
        {
            return ViewHost(points_.view_host(), weights_.view_host());
        }


    private:
        DualPointView points_;
        DualWeightView weights_;
    };

}

#endif //UTOPIA_KOKKOS_QUADRATURE_VIEW_HPP
