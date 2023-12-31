#ifndef UTOPIA_ASSEMBLY_VIEW_HPP
#define UTOPIA_ASSEMBLY_VIEW_HPP

#include "utopia_FunctionSpace.hpp"
#include "utopia_ViewForwardDeclarations.hpp"

namespace utopia {

    template <class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class PhysicalPoint {
    public:
        static const int Dim = Elem::Dim;

        UTOPIA_INLINE_FUNCTION PhysicalPoint(const Quadrature &q) : q_(q), elem_(nullptr) {}

        template <class Point>
        UTOPIA_INLINE_FUNCTION void get(const int qp_idx, Point &p) const {
            Point temp;
            q_.point(qp_idx, temp);
            elem_->point(temp, p);
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

    template <class Mesh, int NComponents, class Quadrature, class... Args>
    class PhysicalPoint<FunctionSpace<Mesh, NComponents, Args...>, Quadrature> {
    public:
        using FunctionSpace = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Elem = typename FunctionSpace::ViewDevice::Elem;

        using ViewDevice = utopia::PhysicalPoint<Elem, typename Quadrature::ViewDevice>;

        PhysicalPoint(const FunctionSpace &, const Quadrature &q) : q_(q) {}

        ViewDevice view_device() const { return ViewDevice(q_.view_device()); }

    private:
        const Quadrature &q_;
    };

    template <class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class PhysicalGradient {
    public:
        static const int Dim = Elem::Dim;
        using Scalar = typename Elem::Scalar;
        //        using SizeType = typename Elem::SizeType;
        using Point = typename Elem::Point;
        using GradValue = typename Elem::GradValue;

        UTOPIA_INLINE_FUNCTION PhysicalGradient(const Quadrature &q) : q_(q), elem_(nullptr) {}

        template <class Grad>
        UTOPIA_INLINE_FUNCTION void get(const int fun_num, const int qp_idx, Grad &g) const {
            Point temp;
            q_.point(qp_idx, temp);
            elem_->grad(fun_num, temp, g);
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

    template <class Mesh, int NComponents, class Quadrature, typename... Args>
    class PhysicalGradient<FunctionSpace<Mesh, NComponents, Args...>, Quadrature> {
    public:
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

    //////////////////////////////////////////////////////////////////////////////////

    template <class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class ShapeFunction {
    public:
        static const int Dim = Elem::Dim;
        using Scalar = typename Elem::Scalar;
        using Point = typename Elem::Point;
        using FunValue = typename Elem::FunValue;

        UTOPIA_INLINE_FUNCTION ShapeFunction(const Quadrature &q) : q_(q), elem_(nullptr) {}

        UTOPIA_INLINE_FUNCTION FunValue get(const int fun_num, const int qp_idx) const {
            Point temp;
            q_.point(qp_idx, temp);
            return elem_->fun(fun_num, temp);
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

    template <class Mesh, int NComponents, class Quadrature, typename... Args>
    class ShapeFunction<FunctionSpace<Mesh, NComponents, Args...>, Quadrature> {
    public:
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

    //////////////////////////////////////////////////////////////////////////////////

    template <class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class Differential {
    public:
        static const int Dim = Elem::Dim;
        using Scalar = typename Elem::Scalar;
        using Point = utopia::StaticVector<Scalar, Dim>;

        UTOPIA_INLINE_FUNCTION Differential(const Quadrature &q) : q_(q), elem_(nullptr) {}

        UTOPIA_INLINE_FUNCTION Scalar operator()(const int qp_idx) const { return get(qp_idx); }

        UTOPIA_INLINE_FUNCTION Scalar get(const int qp_idx) const {
            UTOPIA_DEVICE_ASSERT(elem_->is_affine());
            return q_.weight(qp_idx) * elem_->measure();
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

    template <class Mesh, int NComponents, class Quadrature, typename... Args>
    class Differential<FunctionSpace<Mesh, NComponents, Args...>, Quadrature> {
    public:
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

}  // namespace utopia

#endif  // UTOPIA_ASSEMBLY_VIEW_HPP
