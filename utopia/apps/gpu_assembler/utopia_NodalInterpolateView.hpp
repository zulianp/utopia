#ifndef UTOPIA_NODAL_INTERPOLATE_VIEW_HPP
#define UTOPIA_NODAL_INTERPOLATE_VIEW_HPP

#include "utopia_AssemblyView.hpp"
#include "utopia_DeviceView.hpp"

namespace utopia {

    template<
        class FunctionSpaceView,
        class VectorView,
        class ShapeFunView,
        class QuadratureView,
        class Elem = typename FunctionSpaceView::Elem,
        class MemType = typename Elem::MemType, typename...>
    class NodalInterpolateView {
    public:
        static const int Dim = Elem::Dim;
        using Scalar = typename Elem::Scalar;
        using Point  = utopia::StaticVector<Scalar, Dim>;
        using Eval   = utopia::StaticVector<Scalar, QuadratureView::NPoints>;
        using Coeff  = utopia::StaticVector<Scalar, Elem::NNodes>;
        using DofIndex = typename FunctionSpaceView::DofIndex;

        UTOPIA_INLINE_FUNCTION NodalInterpolateView(const FunctionSpaceView &space, const VectorView &vec, const ShapeFunView &fun)
        : space_(space), vec_(vec), fun_(fun), elem_(nullptr)
        {}

        UTOPIA_INLINE_FUNCTION std::size_t size() const
        {
            return fun_.n_points();
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
            Coeff coeff;
            space_.coefficients(elem, vec_, coeff);

            auto shape_i = fun_.make(i, elem);

            const auto n = shape_i.n_points();
            for(SizeType k = 0; k < n; ++k) {
                values[k] = 0.0;
                for(SizeType j = 0; j < shape_i.n_functions(); ++j) {
                    auto f_interp = shape_i(j, k) * coeff(j);
                    values[k] += f_interp * coeff(j);

                }
            }
        }

        const ShapeFunView &fun() const
        {
            return fun_;
        }

    private:
        FunctionSpaceView space_;
        VectorView vec_;
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

        using ViewDevice = utopia::NodalInterpolateView<
                                        typename FunctionSpace::ViewDevice,
                                        DeviceView<const Vector, 1>,
                                        typename ShapeFunction::ViewDevice,
                                        typename Quadrature::ViewDevice
                                        >;

        NodalInterpolate(
            const FunctionSpace &space,
            const Vector &vec,
            const Quadrature &q) : space_(space), vec_(make_ref(vec)), shape_fun_(space_, q) {}

        NodalInterpolate(
            const FunctionSpace &space,
            const Quadrature &q) : space_(space), vec_(nullptr), shape_fun_(space_, q) {}

        ViewDevice view_device() const
        {
            return ViewDevice(space_.view_device(), space_.assembly_view_device(*vec_), shape_fun_.view_device());
        }

        void update(const Vector &vec)
        {
            vec_ = make_ref(vec);
        }

    private:
        const FunctionSpace &space_;
        std::shared_ptr<const Vector> vec_;
        ShapeFunction shape_fun_;
    };

}

#endif //UTOPIA_NODAL_INTERPOLATE_VIEW_HPP
