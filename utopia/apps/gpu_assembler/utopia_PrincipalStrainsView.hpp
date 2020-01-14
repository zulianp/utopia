#ifndef UTOPIA_PRINCIPAL_STRAINS_VIEW_HPP
#define UTOPIA_PRINCIPAL_STRAINS_VIEW_HPP

#include "utopia_LaplacianView.hpp"
#include "utopia_Utils.hpp"
#include "utopia_GradInterpolate.hpp"

namespace utopia {

    template<
        class FunctionSpaceView,
        class GradInterpolateView,
        class Elem = typename FunctionSpaceView::Elem,
        class MemType = typename Elem::MemType,
        typename...
    >
    class PrincipalStrainsView {
    public:
        using Scalar   = typename Elem::Scalar;
        using SizeType = typename Elem::SizeType;

        PrincipalStrainsView(const GradInterpolateView &grad)
        : grad_(grad)
        {}

        template<class Values>
        UTOPIA_INLINE_FUNCTION void get(const std::size_t &qp, const Elem &elem, Values &values) const
        {
            // grad_
        }

        GradInterpolateView grad_;
    };

    template<class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class PrincipalStrains {};

    template<class Mesh, int NComponents, class Quadrature, typename...Args>
    class PrincipalStrains< FunctionSpace<Mesh, NComponents, Args...>, Quadrature> {
    public:
        using FunctionSpace           = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Vector                  = typename FunctionSpace::Vector;
        using GradInterpolate         = utopia::GradInterpolate<FunctionSpace, Quadrature>;

        using FunctionSpaceViewDevice   = typename FunctionSpace::ViewDevice;
        using GradInterpolateViewDevice = typename GradInterpolate::ViewDevice;

        using ViewDevice              = utopia::PrincipalStrainsView<FunctionSpaceViewDevice, GradInterpolateViewDevice>;


        PrincipalStrains(const FunctionSpace &space, const Quadrature &q)
        : grad_(space, q)
        {}

        inline ViewDevice view_device() const
        {
            return ViewDevice(
                grad_.view_device()
            );
        }

        inline void update(const Vector &x)
        {
            vec_ = utopia::make_ref(x);
        }

    private:
        GradInterpolate grad_;
        std::shared_ptr<const Vector> vec_;
    };

}


#endif //UTOPIA_PRINCIPAL_STRAINS_VIEW_HPP
