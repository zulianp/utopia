#ifndef UTOPIA_SIMD_ASSEMBLY_VIEW_HPP
#define UTOPIA_SIMD_ASSEMBLY_VIEW_HPP

#include "utopia_AssemblyView.hpp"
#include "utopia_LaplacianView.hpp"
#include "utopia_simd_Quadrature.hpp"

namespace utopia {

    ////////////////////////////////////////////////////////////////////////////////////////

    template <class Elem, typename T, int Dim_>
    class PhysicalGradient<Elem, simd::Quadrature<T, Dim_>, typename Elem::MemType> {
    public:
        static const int Dim = Dim_;
        using Scalar = typename Elem::Scalar;
        using GradValue = typename Elem::GradValue;
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
        Quadrature q_;
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
        Quadrature q_;
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
        Quadrature q_;
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
                    mat(j, j) += dot(g_test, g_test) * dx(k);

                    for (int l = j + 1; l < n_fun; ++l) {
                        const auto v = dot(g_test, grad(l, k)) * dx(k);
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

}  // namespace utopia

#endif  // UTOPIA_SIMD_ASSEMBLY_VIEW_HPP
