#ifndef UTOPIA_NODAL_INTERPOLATE_VIEW_HPP
#define UTOPIA_NODAL_INTERPOLATE_VIEW_HPP

#include "utopia_AssemblyView.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_Coefficient.hpp"

namespace utopia {

    template<
        class FunctionSpaceView,
        class CoefficientView,
        class ShapeFunView,
        class QuadratureView,
        class Elem = typename FunctionSpaceView::Elem,
        class MemType = typename Elem::MemType, typename...>
    class NodalInterpolateView {
    public:
        static const int Dim = Elem::Dim;
        using Scalar   = typename Elem::Scalar;
        using SizeType = typename Elem::SizeType;
        using Point    = utopia::StaticVector<Scalar, Dim>;
        using Eval     = utopia::StaticVector<typename Elem::FunValue, QuadratureView::NPoints>;
        using Coeff    = utopia::StaticVector<Scalar, Elem::NFunctions>;

        UTOPIA_INLINE_FUNCTION NodalInterpolateView(const CoefficientView &coeff, const ShapeFunView &fun)
        : coeff_(coeff), fun_(fun), elem_(nullptr)
        {}

        UTOPIA_INLINE_FUNCTION std::size_t size() const
        {
            return fun_.n_points();
        }

        UTOPIA_INLINE_FUNCTION Eval make(const Elem &elem) const
        {
            Eval values;
            get(elem, values);
            return values;
        }

        template<class Values>
        UTOPIA_INLINE_FUNCTION void get(const Elem &elem, Values &values) const
        {
            Coeff elem_coeff;
            coeff_.get(elem, elem_coeff);

            auto shape_i = fun_.make(elem);

            const auto n = shape_i.n_points();
            for(SizeType k = 0; k < n; ++k) {
                values[k] = shape_i(0, k) * elem_coeff(0);

                for(SizeType j = 1; j < shape_i.n_functions(); ++j) {
                    values[k] += shape_i(j, k) * elem_coeff(j);
                }
            }
        }

        const ShapeFunView &fun() const
        {
            return fun_;
        }

    private:
        CoefficientView coeff_;
        ShapeFunView fun_;
        const Elem *elem_;
    };

    template<class FunctionSpace, class Quadrature, typename...Args>
    class NodalInterpolate {};

    template<class Mesh, int NComponents, class Quadrature, typename...Args>
    class NodalInterpolate< FunctionSpace<Mesh, NComponents, Args...>, Quadrature> {
    public:
        using FunctionSpace = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Elem = typename FunctionSpace::ViewDevice::Elem;
        using ShapeFunction = utopia::ShapeFunction<FunctionSpace, Quadrature>;

        using Vector     = typename FunctionSpace::Vector;
        using CoefficientViewDevice = typename utopia::Coefficient<FunctionSpace>::ViewDevice;
        using FunctionSpaceViewDevice = typename FunctionSpace::ViewDevice;

        using ViewDevice = utopia::NodalInterpolateView<
                                        FunctionSpaceViewDevice,
                                        CoefficientViewDevice,
                                        typename ShapeFunction::ViewDevice,
                                        typename Quadrature::ViewDevice
                                        >;
        NodalInterpolate(
            const FunctionSpace &space,
            const Quadrature &q) : coeff_(space), shape_fun_(space, q) {}

        ViewDevice view_device() const
        {
            return ViewDevice(coeff_.view_device(), shape_fun_.view_device());
        }

        void update(const Vector &vec)
        {
            coeff_.update(vec);
        }

    private:
        Coefficient<FunctionSpace> coeff_;
        ShapeFunction shape_fun_;
    };

}

#endif //UTOPIA_NODAL_INTERPOLATE_VIEW_HPP
