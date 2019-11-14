#ifndef UTOPIA_KOKKOS_MATRIX_VIEW_HPP
#define UTOPIA_KOKKOS_MATRIX_VIEW_HPP

#include "utopia_ViewForwardDeclarations.hpp"
#include "utopia_Tensor.hpp"
#include "utopia_Algorithms.hpp"
#include "utopia_DeviceExpression.hpp"

#include <utility>

namespace utopia {

    template<class ArrayView2D_>
    class TensorView<ArrayView2D_, 2> final : public DeviceExpression< TensorView<ArrayView2D_, 2> > {
    public:
        using ArrayView2D = ArrayView2D_;
        using Scalar   = typename Traits<ArrayView2D>::Scalar;
        using SizeType = typename Traits<ArrayView2D>::SizeType;

        enum {
            StoreAs = UTOPIA_BY_REFERENCE
        };

        inline std::string get_class() const override
        {
            return "MatrixView";
        }

        UTOPIA_FUNCTION TensorView() {}
        UTOPIA_FUNCTION TensorView(ArrayView2D &&view) : view_(std::move(view)) {}
        UTOPIA_FUNCTION TensorView(const ArrayView2D &view) : view_(view) {}

        UTOPIA_FUNCTION TensorView(TensorView &&other) : view_(std::move(other.view_)) {}

        template<class... Args>
        UTOPIA_FUNCTION TensorView(const DelegateArgs &, Args && ...args)
        : view_(std::forward<Args>(args)...)
        {}

        // template<class... Args>
        // UTOPIA_FUNCTION TensorView(Args && ...args)
        // : view_(std::forward<Args>(args)...)
        // {}

        template<class Expr>
        UTOPIA_FUNCTION TensorView(const DeviceExpression<Expr> &expr)
        {
            DeviceAssign<TensorView, Expr>::apply(*this, expr.derived());
        }

        template<class Expr>
        UTOPIA_FUNCTION TensorView(DeviceExpression<Expr> &&expr)
        {
            DeviceAssign<TensorView, Expr>::apply(*this, std::move(expr.derived()));
        }

        template<class Expr>
        UTOPIA_INLINE_FUNCTION TensorView &operator=(const DeviceExpression<Expr> &expr)
        {
            DeviceAssign<TensorView, Expr>::apply(*this, expr.derived());
            return *this;
        }

        template<class Expr>
        UTOPIA_INLINE_FUNCTION TensorView &operator+=(const DeviceExpression<Expr> &expr)
        {
            DeviceInPlace<TensorView, Expr, Plus, 2>::apply(*this, expr.derived());
            return *this;
        }

        template<class Expr>
        UTOPIA_INLINE_FUNCTION TensorView &operator*=(const DeviceExpression<Expr> &expr)
        {
            DeviceInPlace<TensorView, Expr, Multiplies, 2>::apply(*this, expr.derived());
            return *this;
        }

        template<class Expr>
        UTOPIA_INLINE_FUNCTION TensorView &operator/=(const DeviceExpression<Expr> &expr)
        {
            DeviceInPlace<TensorView, Expr, Divides, 2>::apply(*this, expr.derived());
            return *this;
        }

        /////////////////////////////////////////////////////////////

        UTOPIA_INLINE_FUNCTION TensorView &operator+=(const Scalar &expr)
        {
            device::shift(expr, view_);
            return *this;
        }

        UTOPIA_INLINE_FUNCTION TensorView &operator*=(const Scalar &expr)
        {
            scale(expr);
            return *this;
        }

        UTOPIA_INLINE_FUNCTION TensorView &operator/=(const Scalar &expr)
        {
            scale(1./expr);
            return *this;
        }

        /////////////////////////////////////////////////////////////

        template<class OtherArrayView2D>
        UTOPIA_FUNCTION void copy(const TensorView<OtherArrayView2D, 2> &other)
        {
            UTOPIA_DEVICE_ASSERT(rows() == other.rows());
            UTOPIA_DEVICE_ASSERT(cols() == other.cols());

            device::copy(other.view_, view_);
        }

        UTOPIA_INLINE_FUNCTION ArrayView2D &raw_type() { return view_; }
        UTOPIA_INLINE_FUNCTION const ArrayView2D &raw_type() const { return view_; }

        UTOPIA_INLINE_FUNCTION SizeType rows() const
        {
            return device::extent(view_, 0);
        }

        UTOPIA_INLINE_FUNCTION SizeType cols() const
        {
            return device::extent(view_, 1);
        }

        UTOPIA_INLINE_FUNCTION SizeType size() const
        {
            return view_.size();
        }

        UTOPIA_INLINE_FUNCTION Scalar &operator()(const SizeType &i, const SizeType &j)
        {
            return view_(i, j);
        }

        UTOPIA_INLINE_FUNCTION const Scalar &operator()(const SizeType &i, const SizeType &j) const
        {
            return view_(i, j);
        }

        UTOPIA_INLINE_FUNCTION const Scalar &get(const SizeType &i, const SizeType &j) const
        {
            return view_(i, j);
        }

        UTOPIA_INLINE_FUNCTION void set(const SizeType &i, const SizeType &j, const Scalar &value)
        {
            view_(i, j) = value;
        }

        UTOPIA_INLINE_FUNCTION void add(const SizeType &i, const SizeType &j, const Scalar &value)
        {
            view_(i, j) += value;
        }

        UTOPIA_INLINE_FUNCTION void scale(const Scalar &alpha)
        {
            device::scale(alpha, view_);
        }

        UTOPIA_INLINE_FUNCTION void axpy(const Scalar &alpha, const TensorView &x)
        {
            device::axpy(alpha, x.view_, view_);
        }

        UTOPIA_INLINE_FUNCTION Scalar dot(const TensorView &other) const
        {
            return device::dot(view_, other.view_);
        }

        template<class RightView, class ResultView>
        UTOPIA_INLINE_FUNCTION void multiply(
            const VectorView<RightView> &right,
            VectorView<ResultView> &result) const
        {
            UTOPIA_DEVICE_ASSERT(!result.is_alias(right));
            UTOPIA_DEVICE_ASSERT(result.size() == rows());
            UTOPIA_DEVICE_ASSERT(right.size() == cols());

            device::mv(view_, right.raw_type(), result.raw_type());
        }

        template<class RightView, class ResultView>
        UTOPIA_INLINE_FUNCTION void multiply(const TensorView<RightView, 2> &right, TensorView<ResultView, 2> &result) const
        {
            UTOPIA_DEVICE_ASSERT(!result.is_alias(right));
            UTOPIA_DEVICE_ASSERT(!is_alias(result));
            UTOPIA_DEVICE_ASSERT(rows() == result.rows());
            UTOPIA_DEVICE_ASSERT(right.cols() == result.cols());
            UTOPIA_DEVICE_ASSERT(cols() == right.rows());

            device::mm(view_, right.view_, result.view_);
        }

        UTOPIA_INLINE_FUNCTION void set(const Scalar &alpha)
        {
            device::fill(alpha, view_);
        }

        UTOPIA_INLINE_FUNCTION void wrap(const ArrayView2D &view)
        {
            view_ = view;
        }

        UTOPIA_INLINE_FUNCTION bool is_alias(const TensorView &other) const
        {
            return &(view_(0, 0)) == &(other.view_(0, 0));
        }

        template<class OtherArrayView2D>
        UTOPIA_INLINE_FUNCTION constexpr static bool is_alias(const TensorView<OtherArrayView2D, 2> &)
        {
            return false;
        }

        inline void describe() const
        {
            const SizeType r = rows();
            const SizeType c = cols();

            for(SizeType i = 0; i < r; ++i) {
                for(SizeType j = 0; j < c; ++j) {
                    std::cout << get(i, j) << "\t";
                }

                std::cout << "\n";
            }
        }

    private:
        ArrayView2D view_;

        UTOPIA_FUNCTION TensorView(const TensorView &) {
            UTOPIA_DEVICE_ASSERT(false);
        }
    };

}

#endif //UTOPIA_KOKKOS_MATRIX_VIEW_HPP
