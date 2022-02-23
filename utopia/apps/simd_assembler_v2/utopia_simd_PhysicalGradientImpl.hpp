#ifndef UTOPIA_SIMD_PHYSICALGRADIENT_IMPL_HPP
#define UTOPIA_SIMD_PHYSICALGRADIENT_IMPL_HPP

#include "utopia_AssemblyView.hpp"

namespace utopia {

    template <class Elem, class Quadrature, class MemType = typename Elem::MemType>
    class SIMDPhysicalGradientImpl_Elem {
    public:
        static const int Dim = Elem::Dim;
        using Scalar = typename simd_v2::FETraits<Elem>::Scalar;
        using GradValue = typename simd_v2::FETraits<Elem>::GradValue;

        UTOPIA_INLINE_FUNCTION SIMDPhysicalGradientImpl_Elem(const Quadrature &q) : q_(q), elem_(nullptr) {}

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

        UTOPIA_INLINE_FUNCTION SIMDPhysicalGradientImpl_Elem make(const Elem &elem) const {
            SIMDPhysicalGradientImpl_Elem pp(q_);
            pp.elem_ = &elem;
            return pp;
        }

    private:
        const Quadrature &q_;
        const Elem *elem_;
    };

    template <class FunctionSpace, class Quadrature, class MemType = typename FunctionSpace::Elem::MemType>
    class SIMDPhysicalGradientImpl_Space {
    public:
        using Elem = typename FunctionSpace::ViewDevice::Elem;
        using ViewDevice = utopia::SIMDPhysicalGradientImpl_Elem<Elem, typename Quadrature::ViewDevice, MemType>;
        using ViewHost = ViewDevice;

        SIMDPhysicalGradientImpl_Space(const FunctionSpace &, const Quadrature &q) : q_(q) {}

        ViewDevice view_device() const { return ViewDevice(q_.view_device()); }

        ViewHost view_host() const { return ViewHost(q_.view_host()); }

    private:
        const Quadrature &q_;
    };

    ////////////////////////////////////////////////////////

    // Uniform gradient
    template <class Elem, class Quadrature>
    class SIMDPhysicalGradientImpl_Elem<Elem, Quadrature, Uniform<>> {
    public:
        static const int Dim = Elem::Dim;
        using Scalar = typename simd_v2::FETraits<Elem>::Scalar;
        using GradValue = typename simd_v2::FETraits<Elem>::GradValue;

        static const int NFunctions = Elem::NFunctions;

        UTOPIA_INLINE_FUNCTION SIMDPhysicalGradientImpl_Elem(const Vc::vector<GradValue> &values) : values_(values) {}

        template <class Grad>
        UTOPIA_INLINE_FUNCTION void get(const int fun_num, const int qp_idx, Grad &g) const {
            g = values_[qp_idx * NFunctions + fun_num];
        }

        UTOPIA_INLINE_FUNCTION const GradValue &operator()(const int fun_num, const int qp_idx) const {
            return values_[qp_idx * NFunctions + fun_num];
        }

        UTOPIA_INLINE_FUNCTION std::size_t n_points() const { return values_.size() / NFunctions; }

        UTOPIA_INLINE_FUNCTION constexpr static std::size_t n_functions() { return NFunctions; }

        UTOPIA_INLINE_FUNCTION SIMDPhysicalGradientImpl_Elem make(const Elem &) const {
            SIMDPhysicalGradientImpl_Elem pp(values_);
            return pp;
        }

    private:
        const Vc::vector<GradValue> &values_;
    };

    template <class FunctionSpace, class Quadrature>
    class SIMDPhysicalGradientImpl_Space<FunctionSpace, Quadrature, Uniform<>> {
    public:
        using Elem = typename FunctionSpace::ViewDevice::Elem;
        static const int NFunctions = Elem::NFunctions;
        using GradValue = typename simd_v2::FETraits<Elem>::GradValue;

        using ViewDevice = utopia::PhysicalGradient<Elem, typename Quadrature::ViewDevice>;
        using ViewHost = utopia::PhysicalGradient<Elem, typename Quadrature::ViewHost>;

        SIMDPhysicalGradientImpl_Space(const FunctionSpace &space, const Quadrature &q) { init(space, q); }

        ViewDevice view_device() const { return ViewDevice(values_); }

        ViewHost view_host() const { return ViewHost(values_); }

    private:
        Vc::vector<GradValue> values_;

        void init(const FunctionSpace &space, const Quadrature &q) {
            Elem e;
            space.elem(0, e);

            const std::size_t n_qp = q.n_points();
            assert(n_qp > 0);

            values_.resize(n_qp * NFunctions);

            for (std::size_t qp = 0; qp < n_qp; ++qp) {
                for (std::size_t i = 0; i < NFunctions; ++i) {
                    e.grad(i, q.point(qp), values_[qp * NFunctions + i]);
                }
            }
        }
    };

}  // namespace utopia

#endif  // UTOPIA_SIMD_PHYSICALGRADIENT_IMPL_HPP
