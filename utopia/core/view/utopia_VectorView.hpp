#ifndef UTOPIA_KOKKOS_VECTOR_VIEW_HPP
#define UTOPIA_KOKKOS_VECTOR_VIEW_HPP

#include "utopia_Accessor.hpp"
#include "utopia_Algorithms.hpp"
#include "utopia_DeviceExpression.hpp"
#include "utopia_TensorView.hpp"
#include "utopia_ViewForwardDeclarations.hpp"
#include "utopia_ViewTraits.hpp"

namespace utopia {

    template <class View>
    class Accessor<TensorView<View, 1>> {
    public:
        using Scalar = typename Traits<TensorView<View, 1>>::Scalar;
        using SizeType = typename Traits<TensorView<View, 1>>::SizeType;

        UTOPIA_INLINE_FUNCTION static const Scalar &get(const TensorView<View, 1> &t, const SizeType &i) {
            return t(i);
        }

        UTOPIA_INLINE_FUNCTION static void set(TensorView<View, 1> &t, const SizeType &i, const Scalar &val) {
            t(i) = val;
        }
    };

    template <class ArrayView_>
    class TensorView<ArrayView_, 1> final : public DeviceExpression<TensorView<ArrayView_, 1>> {
    public:
        using ArrayView = ArrayView_;
        using Scalar = typename Traits<ArrayView>::Scalar;
        using SizeType = typename Traits<ArrayView>::SizeType;

        enum { StoreAs = UTOPIA_BY_REFERENCE };

        inline std::string get_class() const override { return "VectorView"; }

        template <class Expr>
        UTOPIA_INLINE_FUNCTION TensorView &operator=(const DeviceExpression<Expr> &expr) {
            DeviceAssign<TensorView, Expr>::apply(*this, expr.derived());
            return *this;
        }

        /////////////////////////////////////////////////////////////

        template <class Expr>
        UTOPIA_INLINE_FUNCTION TensorView &operator+=(const DeviceExpression<Expr> &expr) {
            DeviceInPlace<TensorView, Expr, Plus, 1>::apply(*this, expr.derived());
            return *this;
        }

        template <class Expr>
        UTOPIA_INLINE_FUNCTION TensorView &operator-=(const DeviceExpression<Expr> &expr) {
            DeviceInPlace<TensorView, Expr, Minus, 1>::apply(*this, expr.derived());
            return *this;
        }

        template <class Expr>
        UTOPIA_INLINE_FUNCTION TensorView &operator*=(const DeviceExpression<Expr> &expr) {
            DeviceInPlace<TensorView, Expr, Multiplies, 1>::apply(*this, expr.derived());
            return *this;
        }

        template <class Expr>
        UTOPIA_INLINE_FUNCTION TensorView &operator/=(const DeviceExpression<Expr> &expr) {
            DeviceInPlace<TensorView, Expr, Divides, 1>::apply(*this, expr.derived());
            return *this;
        }

        /////////////////////////////////////////////////////////////

        UTOPIA_INLINE_FUNCTION TensorView &operator+=(const Scalar &expr) {
            device::shift(expr, view_);
            return *this;
        }

        UTOPIA_INLINE_FUNCTION TensorView &operator*=(const Scalar &expr) {
            scale(expr);
            return *this;
        }

        UTOPIA_INLINE_FUNCTION TensorView &operator/=(const Scalar &expr) {
            scale(1. / expr);
            return *this;
        }

        /////////////////////////////////////////////////////////////

        template <class OtherArrayView>
        UTOPIA_FUNCTION void copy(const TensorView<OtherArrayView, 1> &other) {
            UTOPIA_DEVICE_ASSERT(size() == other.size());
            device::copy(other.view_, view_);
        }

        UTOPIA_INLINE_FUNCTION ArrayView &raw_type() { return view_; }
        UTOPIA_INLINE_FUNCTION constexpr const ArrayView &raw_type() const { return view_; }

        UTOPIA_INLINE_FUNCTION constexpr SizeType size() const { return view_.size(); }

        UTOPIA_INLINE_FUNCTION Scalar &operator()(const SizeType &i) { return view_[i]; }

        UTOPIA_INLINE_FUNCTION constexpr const Scalar &operator()(const SizeType &i) const { return view_[i]; }

        UTOPIA_INLINE_FUNCTION Scalar &operator[](const SizeType &i) { return view_[i]; }

        UTOPIA_INLINE_FUNCTION constexpr const Scalar &operator[](const SizeType &i) const { return view_[i]; }

        UTOPIA_INLINE_FUNCTION constexpr const Scalar &get(const SizeType &i) const { return view_[i]; }

        UTOPIA_INLINE_FUNCTION void set(const SizeType &i, const Scalar &value) { view_[i] = value; }

        UTOPIA_INLINE_FUNCTION void add(const SizeType &i, const Scalar &value) { view_[i] += value; }

        UTOPIA_INLINE_FUNCTION void scale(const Scalar &alpha) { device::scale(alpha, view_); }

        template <class OtherArrayView>
        UTOPIA_INLINE_FUNCTION void axpy(const Scalar &alpha, const TensorView<OtherArrayView, 1> &x) {
            return device::axpy(alpha, x.view_, view_);
        }

        template <class OtherArrayView>
        UTOPIA_INLINE_FUNCTION Scalar dot(const TensorView<OtherArrayView, 1> &other) const {
            return device::dot(view_, other.view_);
        }

        UTOPIA_INLINE_FUNCTION void set(const Scalar &alpha) { device::fill(alpha, view_); }

        UTOPIA_INLINE_FUNCTION void wrap(const ArrayView &view) { view_ = view; }

        UTOPIA_INLINE_FUNCTION constexpr bool is_alias(const TensorView &other) const {
            return &(view_[0]) == &(other.view_[0]);
        }

        template <class OtherArrayView>
        UTOPIA_INLINE_FUNCTION constexpr static bool is_alias(const TensorView<OtherArrayView, 1> &) {
            return false;
        }

        template <class OtherArrayView>
        UTOPIA_INLINE_FUNCTION constexpr bool equals(const TensorView<OtherArrayView, 1> &other,
                                                     const Scalar &tol) const {
            if (size() != other.size()) return false;
            return device::approxeq(view_, other.raw_type(), tol);
        }

        inline void describe(std::ostream &os = std::cout) const {
            const SizeType n = size();
            for (SizeType i = 0; i < n; ++i) {
                disp(get(i), os);
            }
        }

        constexpr TensorView() = default;
        UTOPIA_FUNCTION constexpr TensorView(ArrayView &&view) : view_(std::move(view)) {}
        UTOPIA_FUNCTION constexpr TensorView(const ArrayView &view) : view_(view) {}
        UTOPIA_FUNCTION constexpr TensorView(TensorView &&other) : view_(std::move(other.view_)) {}
        UTOPIA_FUNCTION constexpr TensorView(const TensorView &other) : view_(other.view_) {}

        template <class... Args>
        UTOPIA_FUNCTION TensorView(const DelegateArgs &, Args &&... args) : view_(std::forward<Args>(args)...) {}

        template <class Expr>
        UTOPIA_FUNCTION TensorView(const DeviceExpression<Expr> &expr) {
            DeviceAssign<TensorView, Expr>::apply(*this, expr.derived());
        }

        template <class Expr>
        UTOPIA_FUNCTION TensorView(DeviceExpression<Expr> &&expr) {
            DeviceAssign<TensorView, Expr>::apply(*this, std::move(expr.derived()));
        }

        // private:
        ArrayView view_{};
    };

    template <class View>
    UTOPIA_INLINE_FUNCTION typename Traits<View>::SizeType size(const TensorView<View, 1> &t) {
        return t.size();
    }

    template <class View, typename T, size_t Size>
    inline void convert(const VectorView<View> &in, std::array<T, Size> &out) {
        std::copy(&in[0], &in[0] + Size, &out[0]);
    }
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_VECTOR_VIEW_HPP
