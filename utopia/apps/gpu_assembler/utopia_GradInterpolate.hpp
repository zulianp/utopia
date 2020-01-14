#ifndef UTOPIA_GRAD_INTERPOLATE_HPP
#define UTOPIA_GRAD_INTERPOLATE_HPP

#include "utopia_AssemblyView.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_Coefficient.hpp"

namespace utopia {

    template<
        class FunctionSpaceView,
        class CoefficientView,
        class GradView,
        class QuadratureView,
        class Elem = typename FunctionSpaceView::Elem,
        class MemType = typename Elem::MemType, typename...>
    class GradInterpolateView {
    public:
        static const int Dim = Elem::Dim;
        using Scalar   = typename Elem::Scalar;
        using SizeType = typename Elem::SizeType;
        using Point    = utopia::StaticVector<Scalar, Dim>;
        using Eval     = utopia::StaticVector<typename Elem::GradValue, QuadratureView::NPoints>;
        using Coeff    = utopia::StaticVector<Scalar, Elem::NFunctions>;

        UTOPIA_INLINE_FUNCTION GradInterpolateView(const CoefficientView &coeff, const GradView &grad)
        : coeff_(coeff), grad_(grad), elem_(nullptr)
        {}

        UTOPIA_INLINE_FUNCTION std::size_t size() const
        {
            return grad_.n_points();
        }

        UTOPIA_INLINE_FUNCTION Eval make(const std::size_t &i, const Elem &elem) const
        {
            Eval values;
            get(i, elem, values);
            return values;
        }

        template<class Values>
        UTOPIA_INLINE_FUNCTION void get(const std::size_t &i, const Elem &elem, Values &values) const
        {
            Coeff elem_coeff;
            coeff_.get(elem, elem_coeff);

            auto grad_i = grad_.make(i, elem);

            const auto n = grad_i.n_points();
            for(SizeType k = 0; k < n; ++k) {
                values[k] = grad_i(0, k) * elem_coeff(0);

                for(SizeType j = 1; j < grad_i.n_functions(); ++j) {
                    values[k] += grad_i(j, k) * elem_coeff(j);
                }
            }
        }

        const GradView &grad() const
        {
            return grad_;
        }

    private:
        CoefficientView coeff_;
        GradView grad_;
        const Elem *elem_;
    };


    template<class Space, class Quadrature>
    class GradInterpolate {};


    template<class Mesh, int NComponents, class Quadrature, typename...Args>
    class GradInterpolate< FunctionSpace<Mesh, NComponents, Args...>, Quadrature> {
    public:
        using FunctionSpace = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Elem = typename FunctionSpace::ViewDevice::Elem;
        using PhysicalGradient = utopia::PhysicalGradient<FunctionSpace, Quadrature>;

        using Vector     = typename FunctionSpace::Vector;
        using CoefficientViewDevice = typename utopia::Coefficient<FunctionSpace>::ViewDevice;
        using FunctionSpaceViewDevice = typename FunctionSpace::ViewDevice;

        using ViewDevice = utopia::GradInterpolateView<
                                        FunctionSpaceViewDevice,
                                        CoefficientViewDevice,
                                        typename PhysicalGradient::ViewDevice,
                                        typename Quadrature::ViewDevice
                                        >;
        GradInterpolate(
            const FunctionSpace &space,
            const Quadrature &q) : coeff_(space), grad_(space, q) {}

        ViewDevice view_device() const
        {
            return ViewDevice(coeff_.view_device(), grad_.view_device());
        }

        void update(const Vector &vec)
        {
            coeff_.update(vec);
        }

    private:
        Coefficient<FunctionSpace> coeff_;
        PhysicalGradient grad_;
    };
}

#endif //UTOPIA_GRAD_INTERPOLATE_HPP
