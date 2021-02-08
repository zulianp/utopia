#ifndef UTOPIA_SIMD_ASSEMBLY_VIEW_HPP
#define UTOPIA_SIMD_ASSEMBLY_VIEW_HPP

#include "utopia_AssemblyView.hpp"
#include "utopia_CoefStrainView.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_LaplacianView.hpp"
#include "utopia_NodalInterpolateView.hpp"
#include "utopia_PrincipalShapeStressView.hpp"
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

    template <class FunctionSpaceView,
              class CoefficientView,
              class ShapeFunView,
              class QuadratureView,
              class Elem = typename FunctionSpaceView::Elem,
              class MemType = typename Elem::MemType,
              typename...>
    class VcNodalInterpolateView {
    public:
        static const int Dim = Elem::Dim;
        using Scalar = typename Elem::Scalar;
        // using SizeType = typename Elem::SizeType;
        using Point = utopia::StaticVector<Scalar, Dim>;
        // using Eval = utopia::StaticVector<typename Elem::FunValue, QuadratureView::NPoints>;
        using Coeff = utopia::StaticVector<Scalar, Elem::NFunctions>;

        UTOPIA_INLINE_FUNCTION VcNodalInterpolateView(const CoefficientView &coeff, const ShapeFunView &fun)
            : coeff_(coeff), fun_(fun), elem_(nullptr) {}

        UTOPIA_INLINE_FUNCTION std::size_t size() const { return fun_.n_points(); }

        template <class Values>
        UTOPIA_INLINE_FUNCTION void get(const Elem &elem, Values &values) const {
            Coeff elem_coeff;
            coeff_.get(elem, elem_coeff);

            auto &&shape_i = fun_.make(elem);

            const std::size_t n = shape_i.n_points();
            for (std::size_t k = 0; k < n; ++k) {
                values[k] = shape_i(0, k) * elem_coeff(0);

                for (std::size_t j = 1; j < shape_i.n_functions(); ++j) {
                    values[k] += shape_i(j, k) * elem_coeff(j);
                }
            }
        }

        const ShapeFunView &fun() const { return fun_; }

    private:
        CoefficientView coeff_;
        ShapeFunView fun_;
        const Elem *elem_;
    };

    template <class Mesh, int NComponents, typename T, int Dim_, typename... Args>
    class NodalInterpolate<FunctionSpace<Mesh, NComponents, Args...>, simd::Quadrature<T, Dim_>> {
    public:
        using Quadrature = simd::Quadrature<T, Dim_>;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Elem = typename FunctionSpace::ViewDevice::Elem;
        using ShapeFunction = utopia::ShapeFunction<FunctionSpace, Quadrature>;

        using Vector = typename FunctionSpace::Vector;
        using CoefficientViewDevice = typename utopia::Coefficient<FunctionSpace>::ViewDevice;
        using FunctionSpaceViewDevice = typename FunctionSpace::ViewDevice;
        using Coefficient = utopia::Coefficient<FunctionSpace>;

        using ViewDevice = utopia::VcNodalInterpolateView<FunctionSpaceViewDevice,
                                                          CoefficientViewDevice,
                                                          typename ShapeFunction::ViewDevice,
                                                          typename Quadrature::ViewDevice>;
        NodalInterpolate(const FunctionSpace &space, const Quadrature &q)
            : coeff_(std::make_shared<Coefficient>(space)), shape_fun_(space, q) {}

        NodalInterpolate(const std::shared_ptr<Coefficient> &coeff, const Quadrature &q)
            : coeff_(coeff), shape_fun_(coeff->space(), q) {}

        ViewDevice view_device() const { return ViewDevice(coeff_->view_device(), shape_fun_.view_device()); }

        void update(const Vector &vec) { coeff_->update(vec); }

    private:
        std::shared_ptr<Coefficient> coeff_;
        ShapeFunction shape_fun_;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    template <class FunctionSpaceView,
              class CoefficientView,
              class GradView,
              typename T,
              int Dim_,
              class Elem,
              class MemType,
              typename... Args>
    class GradInterpolateView<FunctionSpaceView,
                              CoefficientView,
                              GradView,
                              simd::Quadrature<T, Dim_>,
                              Elem,
                              MemType,
                              Args...> {
    public:
        using Quadrature = simd::Quadrature<T, Dim_>;
        static const int Dim = Dim_;
        using Scalar = typename Elem::Scalar;
        // using SizeType = typename Elem::SizeType;
        using Point = utopia::StaticVector<Scalar, Dim>;
        using Coeff = utopia::StaticVector<Scalar, Elem::NFunctions>;

        UTOPIA_INLINE_FUNCTION GradInterpolateView(const CoefficientView &coeff, const GradView &grad)
            : coeff_(coeff), grad_(grad), elem_(nullptr) {}

        UTOPIA_INLINE_FUNCTION std::size_t size() const { return grad_.n_points(); }

        template <class Values>
        UTOPIA_INLINE_FUNCTION void get(const Elem &elem, Values &values) const {
            Coeff elem_coeff;
            coeff_.get(elem, elem_coeff);

            auto grad_i = grad_.make(elem);

            const int n = grad_i.n_points();
            for (int k = 0; k < n; ++k) {
                values[k] = grad_i(0, k) * elem_coeff(0);

                for (std::size_t j = 1; j < grad_i.n_functions(); ++j) {
                    values[k] += grad_i(j, k) * elem_coeff(j);
                }
            }
        }

        const GradView &grad() const { return grad_; }

    private:
        CoefficientView coeff_;
        GradView grad_;
        const Elem *elem_;
    };

    ////////////////////////////////////////////////////////////////////////////////////////

    template <class FunctionSpaceView,
              class GradInterpolateView,
              class Elem = typename FunctionSpaceView::Elem,
              class MemType = typename Elem::MemType,
              typename...>
    class VcCoefStrainView {
    public:
        static const int Dim = Elem::Dim;

        using Scalar = typename Elem::Scalar;
        using GradValue = typename Elem::GradValue;

        VcCoefStrainView(const GradInterpolateView &grad) : grad_(grad) {}

        template <class Strain>
        UTOPIA_INLINE_FUNCTION void get(const Elem &elem, Strain &strain) const {
            grad_.get(elem, strain);

            const SizeType n_qp = grad_.size();
            for (SizeType qp = 0; qp < n_qp; ++qp) {
                strain[qp].symmetrize();
            }
        }

    private:
        GradInterpolateView grad_;
    };

    template <class Mesh, int NComponents, typename T, int Dim_, typename... Args>
    class CoefStrain<FunctionSpace<Mesh, NComponents, Args...>, simd::Quadrature<T, Dim_>> {
    public:
        using Quadrature = simd::Quadrature<T, Dim_>;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Vector = typename FunctionSpace::Vector;
        using GradInterpolate = utopia::GradInterpolate<FunctionSpace, Quadrature>;

        using FunctionSpaceViewDevice = typename FunctionSpace::ViewDevice;
        using GradInterpolateViewDevice = typename GradInterpolate::ViewDevice;

        using ViewDevice = utopia::VcCoefStrainView<FunctionSpaceViewDevice, GradInterpolateViewDevice>;

        CoefStrain(const std::shared_ptr<Coefficient<FunctionSpace>> &coeff, const Quadrature &q) : grad_(coeff, q) {}

        inline ViewDevice view_device() const { return ViewDevice(grad_.view_device()); }

        inline void update(const Vector &x) { grad_.update(x); }

    private:
        GradInterpolate grad_;
    };

    //////////////////////////////////////////////////////////////////////////////////////////////////////

    template <class Mesh, int NComponents, typename T, int Dim_, typename... Args>
    class ShapeStress<FunctionSpace<Mesh, NComponents, Args...>, simd::Quadrature<T, Dim_>> {
    public:
        using Quadrature = simd::Quadrature<T, Dim_>;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Vector = typename FunctionSpace::Vector;
        using Scalar = typename FunctionSpace::Scalar;

        using FunctionSpaceViewDevice = typename FunctionSpace::ViewDevice;

        using Elem = typename FunctionSpace::ViewDevice::Elem;
        using GradValue = typename simd::FETraits<Elem, T>::GradValue;
        static const int Dim = Elem::Dim;

        static const int NFunctions = Elem::NFunctions;

        class ViewDevice {
        public:
            ViewDevice() = default;
            ViewDevice(ViewDevice &&other) = default;

            ViewDevice(const ViewDevice &other) {
                const int n = other.stress.size();
                stress.resize(n);

                for (int i = 0; i < n; ++i) {
                    stress[i] = other.stress[i];
                }
            }

            inline void resize(const int n_qp) { stress.resize(n_qp * NFunctions); }

            inline GradValue &operator()(const int fun, const int qp) { return stress[qp * NFunctions + fun]; }

            inline const GradValue &operator()(const int fun, const int qp) const {
                return stress[qp * NFunctions + fun];
            }

        private:
            Vc::vector<GradValue> stress;
        };

        ShapeStress(const FunctionSpace &space, const Quadrature &q, const Scalar &mu, const Scalar &lambda)

        {
            compute_aggregate_stress(space, q, mu, lambda);
        }

        inline const ViewDevice &view_device() const { return view_device_; }

    private:
        ViewDevice view_device_;

        void compute_aggregate_stress(const FunctionSpace &space,
                                      const Quadrature &q,
                                      const Scalar &mu,
                                      const Scalar &lambda) {
            PhysicalGradient<FunctionSpace, Quadrature> grad(space, q);
            auto grad_view = grad.view_host();

            int n_qp = q.n_points();
            view_device_.resize(n_qp);

            Elem e;
            space.elem(0, e);

            auto g = grad_view.make(e);

            GradValue strain;

            for (SizeType i = 0; i < e.n_functions(); ++i) {
                for (int k = 0; k < n_qp; ++k) {
                    g.get(i, k, strain);
                    strain.symmetrize();

                    view_device_(i, k) = 2.0 * mu * strain + lambda * trace(strain) * (device::identity<Scalar>());
                }
            }
        }
    };  // namespace utopia

}  // namespace utopia

#endif  // UTOPIA_SIMD_ASSEMBLY_VIEW_HPP
