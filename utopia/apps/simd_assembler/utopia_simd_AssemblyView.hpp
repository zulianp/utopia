#ifndef UTOPIA_SIMD_ASSEMBLY_VIEW_HPP
#define UTOPIA_SIMD_ASSEMBLY_VIEW_HPP

#include "utopia_AssemblyView.hpp"
#include "utopia_LaplacianView.hpp"
#include "utopia_ProjectionView.hpp"
#include "utopia_StrainView.hpp"

#include "utopia_simd_Quadrature.hpp"

namespace utopia {

    namespace simd {
        template <class Elem, typename T>
        struct FETraits {
            static const int Dim = Elem::Dim;
            using Scalar = T;
            using FunValue = Scalar;
            using Point = simd::Vector<T, Dim>;
            using GradValue = simd::Vector<T, Dim>;
            using STGradX = simd::Vector<T, Dim - 1>;
        };

        template <class Elem, int NVar, typename T>
        struct FETraits<MultiVariateElem<Elem, NVar>, T> {
            static const int Dim = Elem::Dim;
            using Scalar = T;
            using FunValue = simd::Vector<T, NVar>;
            using Point = simd::Vector<T, Dim>;
            using GradValue = simd::Matrix<T, NVar, Dim>;
            using STGradX = simd::Vector<T, Dim - 1>;
        };

        template <class Elem, typename T>
        struct FETraits<MultiVariateElem<Elem, 1>, T> : FETraits<Elem, T> {};

    }  // namespace simd

    ////////////////////////////////////////////////////////////////////////////////////////

    template <class Elem, typename T, int Dim_>
    class PhysicalGradient<Elem, simd::Quadrature<T, Dim_>, typename Elem::MemType> {
    public:
        static const int Dim = Dim_;
        using Scalar = typename simd::FETraits<Elem, T>::Scalar;
        using GradValue = typename simd::FETraits<Elem, T>::GradValue;
        using Quadrature = simd::Quadrature<T, Dim_>;

        UTOPIA_INLINE_FUNCTION PhysicalGradient(const Quadrature &q) : q_(q), elem_(nullptr) {}

        template <class Grad>
        UTOPIA_INLINE_FUNCTION void get(const int fun_num, const int qp_idx, Grad &g) const {
            elem_->grad(fun_num, q_.point(qp_idx), g);
        }

        UTOPIA_INLINE_FUNCTION GradValue operator()(const int fun_num, const int qp_idx) const {
            GradValue g;
            get(fun_num, qp_idx, g);
            return g;
        }

        UTOPIA_INLINE_FUNCTION std::size_t n_points() const { return q_.n_points(); }

        UTOPIA_INLINE_FUNCTION constexpr static std::size_t n_functions() { return Elem::NFunctions; }

        UTOPIA_INLINE_FUNCTION PhysicalGradient make(const Elem &elem) const {
            PhysicalGradient pp(q_);
            pp.elem_ = &elem;
            return pp;
        }

    private:
        const Quadrature &q_;
        const Elem *elem_;
    };

    template <class Mesh, int NComponents, typename T, int Dim_, typename... Args>
    class PhysicalGradient<FunctionSpace<Mesh, NComponents, Args...>, simd::Quadrature<T, Dim_>> {
    public:
        using Quadrature = simd::Quadrature<T, Dim_>;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Elem = typename FunctionSpace::ViewDevice::Elem;

        using ViewDevice = utopia::PhysicalGradient<Elem, typename Quadrature::ViewDevice>;
        using ViewHost = utopia::PhysicalGradient<Elem, typename Quadrature::ViewHost>;

        PhysicalGradient(const FunctionSpace &, const Quadrature &q) : q_(q) {}

        ViewDevice view_device() const { return ViewDevice(q_.view_device()); }

        ViewHost view_host() const { return ViewHost(q_.view_host()); }

    private:
        const Quadrature &q_;
    };

    ////////////////////////////////////////////////////////////////////////////////////////

    template <class Elem, typename T, int Dim_>
    class Differential<Elem, simd::Quadrature<T, Dim_>, typename Elem::MemType> {
    public:
        static const int Dim = Dim_;
        using Scalar = T;
        using Quadrature = simd::Quadrature<T, Dim_>;

        UTOPIA_INLINE_FUNCTION Differential(const Quadrature &q) : q_(q), elem_(nullptr) {}

        UTOPIA_INLINE_FUNCTION Scalar operator()(const int qp_idx) const { return get(qp_idx); }

        UTOPIA_INLINE_FUNCTION Scalar get(const int qp_idx) const {
            if (elem_->is_affine()) {
                return q_.weight(qp_idx) * elem_->measure();
            } else {
                UTOPIA_DEVICE_ASSERT(false);
                return -1.0;
            }
        }

        UTOPIA_INLINE_FUNCTION std::size_t size() const { return q_.n_points(); }

        UTOPIA_INLINE_FUNCTION Differential make(const Elem &elem) const {
            Differential pp(q_);
            pp.elem_ = &elem;
            return pp;
        }

    private:
        const Quadrature &q_;
        const Elem *elem_;
    };

    template <class Mesh, int NComponents, typename T, int Dim_, typename... Args>
    class Differential<FunctionSpace<Mesh, NComponents, Args...>, simd::Quadrature<T, Dim_>> {
    public:
        using Quadrature = simd::Quadrature<T, Dim_>;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Elem = typename FunctionSpace::ViewDevice::Elem;

        using ViewDevice = utopia::Differential<Elem, typename Quadrature::ViewDevice>;
        using ViewHost = utopia::Differential<Elem, typename Quadrature::ViewHost>;

        Differential(const FunctionSpace &, const Quadrature &q) : q_(q) {}

        ViewDevice view_device() const { return ViewDevice(q_.view_device()); }

        ViewHost view_host() const { return ViewHost(q_.view_host()); }

    private:
        const Quadrature &q_;
    };

    ////////////////////////////////////////////////////////////////////////////////////////

    template <class Elem, typename T, int Dim_>
    class PhysicalPoint<Elem, simd::Quadrature<T, Dim_>, typename Elem::MemType> {
    public:
        using Quadrature = simd::Quadrature<T, Dim_>;
        static const int Dim = Dim_;

        UTOPIA_INLINE_FUNCTION PhysicalPoint(const Quadrature &q) : q_(q), elem_(nullptr) {}

        template <class Point>
        UTOPIA_INLINE_FUNCTION void get(const int qp_idx, Point &p) const {
            elem_->point(q_.point(qp_idx), p);
        }

        UTOPIA_INLINE_FUNCTION std::size_t size() const { return q_.n_points(); }

        UTOPIA_INLINE_FUNCTION PhysicalPoint make(const Elem &elem) const {
            PhysicalPoint pp(q_);
            pp.elem_ = &elem;
            return pp;
        }

    private:
        const Quadrature &q_;
        const Elem *elem_;
    };

    template <class Mesh, int NComponents, typename T, int Dim_, class... Args>
    class PhysicalPoint<FunctionSpace<Mesh, NComponents, Args...>, simd::Quadrature<T, Dim_>> {
    public:
        using Quadrature = simd::Quadrature<T, Dim_>;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Elem = typename FunctionSpace::ViewDevice::Elem;

        using ViewDevice = utopia::PhysicalPoint<Elem, typename Quadrature::ViewDevice>;

        PhysicalPoint(const FunctionSpace &, const Quadrature &q) : q_(q) {}

        ViewDevice view_device() const { return ViewDevice(q_.view_device()); }

    private:
        const Quadrature &q_;
    };

    ////////////////////////////////////////////////////////////////////////////////////////

    template <class Mesh, int NComponents, typename T, int Dim_, typename... Args>
    class Laplacian<FunctionSpace<Mesh, NComponents, Args...>, simd::Quadrature<T, Dim_>> {
    public:
        using Quadrature = simd::Quadrature<T, Dim_>;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Scalar = typename FunctionSpace::Scalar;
        using SizeType = typename FunctionSpace::SizeType;
        using Elem = typename FunctionSpace::ViewDevice::Elem;
        static const int NDofs = FunctionSpace::NDofs;

        using ViewDevice = AssemblerView<StaticMatrix<Scalar, NDofs, NDofs>>;

        Laplacian(const FunctionSpace &space, const Quadrature &q) { init(space, q); }

        ViewDevice view_device() const { return ViewDevice(mat_); }

        template <class Grad, class DX, class Matrix>
        UTOPIA_INLINE_FUNCTION static void assemble(const Grad &grad, const DX &dx, Matrix &mat) {
            const int n = grad.n_points();
            const int n_fun = grad.n_functions();

            for (int k = 0; k < n; ++k) {
                for (int j = 0; j < n_fun; ++j) {
                    const auto g_test = grad(j, k);
                    mat(j, j) += simd::integrate(dot(g_test, g_test) * dx(k));

                    for (int l = j + 1; l < n_fun; ++l) {
                        const auto v = simd::integrate(dot(g_test, grad(l, k)) * dx(k));
                        mat(j, l) += v;
                        mat(l, j) += v;
                    }
                }
            }
        }

    private:
        StaticMatrix<Scalar, NDofs, NDofs> mat_;

        void init(const FunctionSpace &space, const Quadrature &q) {
            PhysicalGradient<FunctionSpace, Quadrature> grad(space, q);
            Differential<FunctionSpace, Quadrature> differential(space, q);

            auto grad_view = grad.view_host();
            auto dx_view = differential.view_host();

            Elem e;
            space.elem(0, e);

            auto g = grad_view.make(e);
            auto dx = dx_view.make(e);

            mat_.set(0.0);
            assemble(g, dx, mat_);
        }
    };

    ////////////////////////////////////////////////////////////////////////////////////////

    template <class Mesh, int NComponents, typename T, int Dim_, class Function, typename... Args>
    class Projection<FunctionSpace<Mesh, NComponents, Args...>, simd::Quadrature<T, Dim_>, Function> {
    public:
        using Quadrature = simd::Quadrature<T, Dim_>;
        using QPoint = typename Quadrature::Point;
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
                        acc(j) += simd::integrate(fun_(p) * shape(j, k) * dx(k));
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

    ////////////////////////////////////////////////////////////////////////////////////////

    template <class Elem, typename T, int Dim_, class MemType>
    class ShapeFunction<Elem, simd::Quadrature<T, Dim_>, MemType> {
    public:
        using Quadrature = simd::Quadrature<T, Dim_>;
        static const int Dim = Dim_;
        using Scalar = typename Elem::Scalar;
        using Point = typename Elem::Point;
        using FunValue = typename simd::FETraits<Elem, T>::FunValue;

        UTOPIA_INLINE_FUNCTION ShapeFunction(const Quadrature &q) : q_(q), elem_(nullptr) {}

        UTOPIA_INLINE_FUNCTION FunValue get(const int fun_num, const int qp_idx) const {
            return elem_->fun(fun_num, q_.point(qp_idx));
        }

        UTOPIA_INLINE_FUNCTION FunValue operator()(const int fun_num, const int qp_idx) const {
            return get(fun_num, qp_idx);
        }

        UTOPIA_INLINE_FUNCTION std::size_t n_points() const { return q_.n_points(); }

        UTOPIA_INLINE_FUNCTION constexpr static std::size_t n_functions() { return Elem::NFunctions; }

        UTOPIA_INLINE_FUNCTION ShapeFunction make(const Elem &elem) const {
            ShapeFunction pp(q_);
            pp.elem_ = &elem;
            return pp;
        }

    private:
        const Quadrature &q_;
        const Elem *elem_;
    };

    template <class Mesh, int NComponents, typename T, int Dim_, typename... Args>
    class ShapeFunction<FunctionSpace<Mesh, NComponents, Args...>, simd::Quadrature<T, Dim_>> {
    public:
        using Quadrature = simd::Quadrature<T, Dim_>;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Elem = typename FunctionSpace::ViewDevice::Elem;

        using ViewDevice = utopia::ShapeFunction<Elem, typename Quadrature::ViewDevice>;
        using ViewHost = utopia::ShapeFunction<Elem, typename Quadrature::ViewHost>;

        ShapeFunction(const FunctionSpace &, const Quadrature &q) : q_(q) {}

        ViewDevice view_device() const { return ViewDevice(q_.view_device()); }

        ViewHost view_host() const { return ViewHost(q_.view_host()); }

    private:
        const Quadrature &q_;
    };

    ////////////////////////////////////////////////////////////////////////////////////////

    template <class Elem, typename T, int Dim_, class MemType, typename... Args>
    class StrainView<Elem, simd::Quadrature<T, Dim_>, MemType, Args...> {
    public:
        using Quadrature = simd::Quadrature<T, Dim_>;
        static const int Dim = Dim_;
        using GradValue = typename simd::FETraits<Elem, T>::GradValue;

        static const std::size_t NFunctions = Elem::NFunctions;

        UTOPIA_INLINE_FUNCTION StrainView(const Vc::vector<GradValue> &values) : values_(values) {}

        UTOPIA_INLINE_FUNCTION const GradValue &operator()(const int fun_num, const int qp_idx) const {
            return values_[qp_idx * NFunctions + fun_num];
        }

        UTOPIA_INLINE_FUNCTION std::size_t n_points() const { return values_.size() / NFunctions; }

        UTOPIA_INLINE_FUNCTION constexpr static std::size_t n_functions() { return NFunctions; }

        UTOPIA_INLINE_FUNCTION const StrainView &make(const Elem & /*elem*/) const { return *this; }

    private:
        const Vc::vector<GradValue> &values_;
    };

    template <class FunctionSpace, typename T, int Dim_>
    class Strain<FunctionSpace, simd::Quadrature<T, Dim_>> {
    public:
        using Quadrature = simd::Quadrature<T, Dim_>;
        using Elem = typename FunctionSpace::ViewDevice::Elem;

        using ViewDevice = utopia::StrainView<Elem, typename Quadrature::ViewDevice>;
        using ViewHost = utopia::StrainView<Elem, typename Quadrature::ViewHost>;

        using GradValue = typename simd::FETraits<Elem, T>::GradValue;

        static const std::size_t NFunctions = Elem::NFunctions;

        Strain(const FunctionSpace &space, const Quadrature &q) : values_() { init(space, q); }

        ViewDevice view_device() const { return ViewDevice(values_); }

    private:
        // ArrayView<GradValue, NQPoints, NFunctions> values_;
        Vc::vector<GradValue> values_;

        void init(const FunctionSpace &space, const Quadrature &q) {
            PhysicalGradient<FunctionSpace, Quadrature> grad(space, q);

            auto grad_view = grad.view_host();

            Elem e;
            space.elem(0, e);

            auto g = grad_view.make(e);

            const std::size_t n_qp = q.n_points();
            values_.resize(n_qp * NFunctions);

            for (std::size_t qp = 0; qp < n_qp; ++qp) {
                for (std::size_t i = 0; i < g.n_functions(); ++i) {
                    auto g_i = g(i, qp);
                    values_[qp * NFunctions + i] = 0.5 * (g_i + transpose(g_i));
                }
            }
        }
    };

    ////////////////////////////////////////////////////////////////////////////////////////

}  // namespace utopia

#endif  // UTOPIA_SIMD_ASSEMBLY_VIEW_HPP
