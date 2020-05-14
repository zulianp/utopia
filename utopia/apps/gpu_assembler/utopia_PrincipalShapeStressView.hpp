#ifndef UTOPIA_PRINCIPAL_SHAPE_STRESS_VIEW_HPP
#define UTOPIA_PRINCIPAL_SHAPE_STRESS_VIEW_HPP

#include "utopia_LaplacianView.hpp"
#include "utopia_Utils.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_split_matrix.hpp"
#include "utopia_DeviceIdentity.hpp"
#include "utopia_TensorView4.hpp"
#include "utopia_DeviceTensorProduct.hpp"
#include "utopia_DeviceTensorContraction.hpp"

namespace utopia {


    template<class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class PrincipalShapeStress {};

    template<class Mesh, int NComponents, class Quadrature, typename...Args>
    class PrincipalShapeStress< FunctionSpace<Mesh, NComponents, Args...>, Quadrature> {
    public:

        using FunctionSpace           = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Vector                  = typename FunctionSpace::Vector;
        using Scalar                  = typename FunctionSpace::Scalar;

        using FunctionSpaceViewDevice   = typename FunctionSpace::ViewDevice;


        using Elem      = typename FunctionSpace::ViewDevice::Elem;
        using GradValue = typename Elem::GradValue;
        static const int Dim = Elem::Dim;

        static const int NFunctions = Elem::NFunctions;

        class ViewDevice {
        public:
            ArrayView<GradValue, NFunctions, Quadrature::NPoints> stress;
            Tensor4th<Scalar, Dim, Dim, Dim, Dim> C;
            ViewDevice() = default;

            ViewDevice(const ViewDevice &other)
            {
                for(int j = 0; j < NFunctions; ++j) {
                    for(int i = 0; i < Quadrature::NPoints; ++i) {
                        stress(j, i).copy(other.stress(j, i));
                    }
                }

                C.copy(other.C); 
            }
        };

        PrincipalShapeStress(
            const FunctionSpace &space,
            const Quadrature &q,
            const Scalar &mu,
            const Scalar &lambda)

        {
            init(space, q, mu, lambda);
        }

        inline const ViewDevice &view_device() const
        {
            return view_device_;
        }

    private:
        ViewDevice view_device_;

        void compute_aggregate_stress(
            const FunctionSpace &space,
            const Quadrature &q,
            const Scalar &mu,
            const Scalar &lambda)
        {
            PhysicalGradient<FunctionSpace, Quadrature> grad(space, q);
            auto grad_view = grad.view_host();

            Elem e;
            space.elem(0, e);

            auto g  = grad_view.make(e);

            GradValue strain;

            for(SizeType i = 0; i < e.n_functions(); ++i) {
                for(SizeType k = 0; k < Quadrature::NPoints; ++k) {
                    g.get(i, k, strain);
                    strain.symmetrize();

                    view_device_.stress(i, k) = 2.0 * mu * strain + lambda * trace(strain) * (device::identity<Scalar>());
                }
            }
        }

        bool kroneckerDelta(const SizeType & i, const SizeType & j){
            return (i==j) ? 1.0 : 0.0; 
        }        

        void fill_in_isotropic_elast_tensor(const Scalar & mu, const Scalar & lambda){
            for(SizeType i = 0; i < Dim; ++i) {
                for(SizeType j = 0; j < Dim; ++j) {
                    for(SizeType k = 0; k < Dim; ++k) {
                        for(SizeType l = 0; l < Dim; ++l) {
                            Scalar val = lambda * kroneckerDelta(i,j)* kroneckerDelta(k,l); 
                            val += mu * (kroneckerDelta(i,k)* kroneckerDelta(j,l)); 
                            val += mu * (kroneckerDelta(i,l)* kroneckerDelta(j,k));
                            view_device_.C.set(i, j, k, l, val);
                        }
                    }
                }
            }            
        }        


        void init(const FunctionSpace &space, const Quadrature &q, const Scalar &mu, const Scalar &lambda)
        {
            compute_aggregate_stress(space, q, mu, lambda);
            fill_in_isotropic_elast_tensor(mu, lambda); 
        }



    };

}

#endif //UTOPIA_PRINCIPAL_SHAPE_STRESS_VIEW_HPP