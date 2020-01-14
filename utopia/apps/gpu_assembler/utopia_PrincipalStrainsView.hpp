#ifndef UTOPIA_PRINCIPAL_STRAINS_VIEW_HPP
#define UTOPIA_PRINCIPAL_STRAINS_VIEW_HPP

#include "utopia_LaplacianView.hpp"
#include "utopia_Utils.hpp"
#include "utopia_GradInterpolate.hpp"

namespace utopia {



    template<class Space, class Quadrature>
    class PrincipalStrainsView {
    public:

    };

    template<class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class PrincipalStrains {
    public:

    };

    template<class Mesh, int NComponents, class Quadrature, typename...Args>
    class PrincipalStrains< FunctionSpace<Mesh, NComponents, Args...>, Quadrature> {
    public:
        using FunctionSpace           = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Vector                  = typename FunctionSpace::Vector;

        using FunctionSpaceViewDevice = typename FunctionSpace::ViewDevice;
        using QuadratureViewDevice    = typename Quadrature::ViewDevice;

        using ViewDevice              = utopia::PrincipalStrainsView<FunctionSpaceViewDevice, QuadratureViewDevice>;


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
        GradInterpolate<FunctionSpace, Quadrature> grad_;
        std::shared_ptr<const Vector> vec_;
    };

}


#endif //UTOPIA_PRINCIPAL_STRAINS_VIEW_HPP
