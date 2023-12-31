#ifndef UTOPIA_PROJECTION_VIEW_HPP
#define UTOPIA_PROJECTION_VIEW_HPP

#include <utility>

#include "utopia_AssemblyView.hpp"
#include "utopia_LaplacianView.hpp"

namespace utopia {

    template <class Elem, class Quadrature, class Function, class MemType = typename Elem::MemType, typename...>
    class Projection {};

    template <class Mesh, int NComponents, class Quadrature, class Function, typename... Args>
    class Projection<FunctionSpace<Mesh, NComponents, Args...>, Quadrature, Function> {
    public:
        using FunctionSpace = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Vector = typename FunctionSpace::Vector;
        using Scalar = typename FunctionSpace::Scalar;
        using Point = typename FunctionSpace::Point;
        using SizeType = typename FunctionSpace::SizeType;
        using Elem = typename FunctionSpace::ViewDevice::Elem;
        static const int NNodes = Elem::NNodes;

        using Differential = utopia::Differential<FunctionSpace, Quadrature>;
        using ShapeFunction = utopia::ShapeFunction<FunctionSpace, Quadrature>;
        using PhysicalPoint = utopia::PhysicalPoint<FunctionSpace, Quadrature>;

        class ViewDevice {
        public:
            using PhysicalPointView = typename PhysicalPoint::ViewDevice;
            using ShapeFunctionView = typename ShapeFunction::ViewDevice;
            using DifferentialView = typename Differential::ViewDevice;

            template <typename SizeType, class Elem, class Accumulator>
            UTOPIA_INLINE_FUNCTION void assemble(const SizeType & /*i*/, const Elem &e, Accumulator &acc) const {
                auto dx = dx_.make(e);
                auto shape = shape_fun_.make(e);
                auto points = point_.make(e);

                Point p;
                const int n = shape.n_points();
                const int n_fun = shape.n_functions();

                for (int k = 0; k < n; ++k) {
                    points.get(k, p);
                    for (int j = 0; j < n_fun; ++j) {
                        acc(j) += fun_(p) * shape(j, k) * dx(k);
                    }
                }
            }

            ViewDevice(Function fun,
                       const PhysicalPointView &points,
                       const ShapeFunctionView &shape_fun,
                       const DifferentialView &dx)
                : fun_(std::move(fun)), point_(points), shape_fun_(shape_fun), dx_(dx) {}

            Function fun_;
            PhysicalPointView point_;
            ShapeFunctionView shape_fun_;
            DifferentialView dx_;
        };

        Projection(const FunctionSpace &space, const Quadrature &q, Function fun)
            : fun_(std::move(fun)), point_(space, q), shape_fun_(space, q), dx_(space, q) {}

        ViewDevice view_device() const {
            return ViewDevice(fun_, point_.view_device(), shape_fun_.view_device(), dx_.view_device());
        }

    private:
        Function fun_;
        PhysicalPoint point_;
        ShapeFunction shape_fun_;
        Differential dx_;
    };
}  // namespace utopia

#endif  // UTOPIA_PROJECTION_VIEW_HPP
