#ifndef UTOPIA_KOKKOS_TENSOR_VIEW_4_HPP
#define UTOPIA_KOKKOS_TENSOR_VIEW_4_HPP

#include "utopia_Accessor.hpp"
#include "utopia_Algorithms.hpp"
#include "utopia_ArrayView4.hpp"
#include "utopia_Base.hpp"
#include "utopia_DeviceExpression.hpp"
#include "utopia_Tensor.hpp"
#include "utopia_TensorView.hpp"
#include "utopia_ViewForwardDeclarations.hpp"

#include <utility>

namespace utopia {

    template <class ArrayView4D_>
    class TensorView<ArrayView4D_, 4> final : public DeviceExpression<TensorView<ArrayView4D_, 4> > {
    public:
        using ArrayView4D = ArrayView4D_;
        using Scalar = typename Traits<ArrayView4D>::Scalar;
        using SizeType = typename Traits<ArrayView4D>::SizeType;

        enum { StoreAs = UTOPIA_BY_REFERENCE };

        inline std::string get_class() const override { return "TensorView4"; }

        UTOPIA_FUNCTION TensorView() = default;
        UTOPIA_FUNCTION TensorView(ArrayView4D &&view) : view_(std::move(view)) {}
        UTOPIA_FUNCTION TensorView(const ArrayView4D &view) : view_(view) {}

        UTOPIA_FUNCTION TensorView(TensorView &&other) : view_(std::move(other.view_)) {}

        template <class... Args>
        UTOPIA_FUNCTION TensorView(const DelegateArgs &, Args &&... args) : view_(std::forward<Args>(args)...) {}

        // template<class... Args>
        // UTOPIA_FUNCTION TensorView(Args && ...args)
        // : view_(std::forward<Args>(args)...)
        // {}

        template <class Expr>
        UTOPIA_FUNCTION TensorView(const DeviceInverse<Expr> &expr) {
            expr.apply(*this);
        }

        template <class Expr>
        UTOPIA_FUNCTION TensorView(const DeviceExpression<Expr> &expr) {
            DeviceAssign<TensorView, Expr>::apply(*this, expr.derived());
        }

        template <class Expr>
        UTOPIA_FUNCTION TensorView(DeviceExpression<Expr> &&expr) {
            DeviceAssign<TensorView, Expr>::apply(*this, std::move(expr.derived()));
        }

        template <class Expr>
        UTOPIA_INLINE_FUNCTION TensorView &operator=(const DeviceExpression<Expr> &expr) {
            DeviceAssign<TensorView, Expr>::apply(*this, expr.derived());
            return *this;
        }

        template <class Expr>
        UTOPIA_INLINE_FUNCTION TensorView &operator+=(const DeviceExpression<Expr> &expr) {
            DeviceInPlace<TensorView, Expr, Plus, 4>::apply(*this, expr.derived());
            return *this;
        }

        template <class Expr>
        UTOPIA_INLINE_FUNCTION TensorView &operator-=(const DeviceExpression<Expr> &expr) {
            DeviceInPlace<TensorView, Expr, Minus, 4>::apply(*this, expr.derived());
            return *this;
        }

        template <class Expr>
        UTOPIA_INLINE_FUNCTION TensorView &operator*=(const DeviceExpression<Expr> &expr) {
            DeviceInPlace<TensorView, Expr, Multiplies, 4>::apply(*this, expr.derived());
            return *this;
        }

        template <class Expr>
        UTOPIA_INLINE_FUNCTION TensorView &operator/=(const DeviceExpression<Expr> &expr) {
            DeviceInPlace<TensorView, Expr, Divides, 4>::apply(*this, expr.derived());
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

        template <class OtherArrayView4D>
        UTOPIA_FUNCTION void copy(const TensorView<OtherArrayView4D, 4> &other) {
            UTOPIA_DEVICE_ASSERT(size() == other.size());
            device::copy(other.raw_type(), view_);
        }

        UTOPIA_INLINE_FUNCTION ArrayView4D &raw_type() { return view_; }
        UTOPIA_INLINE_FUNCTION const ArrayView4D &raw_type() const { return view_; }

        friend UTOPIA_INLINE_FUNCTION SizeType extent(const TensorView &that, const SizeType dim) {
            return device::extent(that.view_, dim);
        }

        UTOPIA_INLINE_FUNCTION SizeType size() const { return view_.size(); }

        UTOPIA_INLINE_FUNCTION Scalar &operator()(const SizeType &i,
                                                  const SizeType &j,
                                                  const SizeType &k,
                                                  const SizeType &l) {
            return view_(i, j, k, l);
        }

        UTOPIA_INLINE_FUNCTION const Scalar &operator()(const SizeType &i,
                                                        const SizeType &j,
                                                        const SizeType &k,
                                                        const SizeType &l) const {
            return view_(i, j, k, l);
        }

        UTOPIA_INLINE_FUNCTION const Scalar &get(const SizeType &i,
                                                 const SizeType &j,
                                                 const SizeType &k,
                                                 const SizeType &l) const {
            return view_(i, j, k, l);
        }

        UTOPIA_INLINE_FUNCTION void set(const SizeType &i,
                                        const SizeType &j,
                                        const SizeType &k,
                                        const SizeType &l,
                                        const Scalar &value) {
            view_(i, j, k, l) = value;
        }

        UTOPIA_INLINE_FUNCTION void add(const SizeType &i,
                                        const SizeType &j,
                                        const SizeType &k,
                                        const SizeType &l,
                                        const Scalar &value) {
            view_(i, j, k, l) += value;
        }

        // template<class Col>
        // UTOPIA_INLINE_FUNCTION void col(const SizeType &j, Col &c) const
        // {
        //     const SizeType r = rows();
        //     UTOPIA_DEVICE_ASSERT(j < cols());

        //     for(SizeType i = 0; i < r; ++i) {
        //         c[i] = get(i, j);
        //     }
        // }

        // template<class Col>
        // UTOPIA_INLINE_FUNCTION void set_col(const SizeType &j, const Col &c)
        // {
        //     const SizeType r = rows();
        //     UTOPIA_DEVICE_ASSERT(j < cols());

        //     for(SizeType i = 0; i < r; ++i) {
        //         set(i, j, c(i));
        //     }
        // }

        UTOPIA_INLINE_FUNCTION void scale(const Scalar &alpha) { device::scale(alpha, view_); }

        UTOPIA_INLINE_FUNCTION void axpy(const Scalar &alpha, const TensorView &x) {
            device::axpy(alpha, x.view_, view_);
        }

        UTOPIA_INLINE_FUNCTION Scalar dot(const TensorView &other) const { return device::dot(view_, other.view_); }

        // template<class RightView, class ResultView>
        // UTOPIA_INLINE_FUNCTION void multiply(
        //     const VectorView<RightView> &right,
        //     VectorView<ResultView> &result) const
        // {
        //     // UTOPIA_DEVICE_ASSERT(!result.is_alias(right));
        //     // UTOPIA_DEVICE_ASSERT(result.size() == rows());
        //     // UTOPIA_DEVICE_ASSERT(right.size() == cols());

        //     device::mv(view_, right.raw_type(), result.raw_type());
        // }

        // template<class RightView, class ResultView>
        // UTOPIA_INLINE_FUNCTION void multiply(const TensorView<RightView, 4> &right, TensorView<ResultView, 4>
        // &result) const
        // {
        //     // UTOPIA_DEVICE_ASSERT(!result.is_alias(right));
        //     // UTOPIA_DEVICE_ASSERT(!is_alias(result));
        //     // UTOPIA_DEVICE_ASSERT(rows() == result.rows());
        //     // UTOPIA_DEVICE_ASSERT(right.cols() == result.cols());
        //     // UTOPIA_DEVICE_ASSERT(cols() == right.rows());

        //     device::mm(view_, right.view_, result.view_);
        // }

        // UTOPIA_INLINE_FUNCTION void shift_diag(const Scalar &alpha)
        // {
        //     device::shift_diag(alpha, view_);
        // }

        // UTOPIA_INLINE_FUNCTION void symmetrize()
        // {
        //     device::symmetrize(view_);
        // }

        // UTOPIA_INLINE_FUNCTION bool is_diagonal(const Scalar &tol) const
        // {
        //     return device::is_diagonal(view_, tol);
        // }

        // UTOPIA_INLINE_FUNCTION bool has_one_nz_per_col(const Scalar &tol) const
        // {
        //     return device::has_one_nz_per_col(view_, tol);
        // }

        UTOPIA_INLINE_FUNCTION void set(const Scalar &alpha) { device::fill(alpha, view_); }

        UTOPIA_INLINE_FUNCTION void wrap(const ArrayView4D &view) { view_ = view; }

        UTOPIA_INLINE_FUNCTION bool is_alias(const TensorView &other) const {
            return &(view_(0, 0)) == &(other.view_(0, 0));
        }

        template <class OtherArrayView4D>
        UTOPIA_INLINE_FUNCTION constexpr static bool is_alias(const TensorView<OtherArrayView4D, 4> &) {
            return false;
        }

        inline void describe(std::ostream &os = std::cout) const {
            const SizeType N0 = extent(*this, 0);
            const SizeType N1 = extent(*this, 1);
            const SizeType N2 = extent(*this, 2);
            const SizeType N3 = extent(*this, 3);

            for (SizeType i = 0; i < N0; ++i) {
                for (SizeType j = 0; j < N1; ++j) {
                    for (SizeType k = 0; k < N2; ++k) {
                        for (SizeType l = 0; l < N3; ++l) {
                            os << "(" << i << "," << j << "," << k << "," << l << ") -> ";
                            os << get(i, j, k, l);
                            os << "\n";
                        }
                    }
                }
            }
        }

        inline void identity(const Scalar &diag_val = 1.0) {
            // https://www.sciencedirect.com/topics/engineering/identity-tensor
            set(0.0);

            const SizeType N0 = extent(*this, 0);
            const SizeType N1 = extent(*this, 1);
            const SizeType N2 = extent(*this, 2);
            const SizeType N3 = extent(*this, 3);

            for (SizeType i = 0; i < N0; ++i) {
                for (SizeType j = 0; j < N1; ++j) {
                    for (SizeType k = 0; k < N2; ++k) {
                        for (SizeType l = 0; l < N3; ++l) {
                            set(i, j, k, l, (i == l) * (j == k) * diag_val);
                        }
                    }
                }
            }
        }

        inline void identity_sym() {
            set(0.0);
            for (SizeType i = 0; i < extent(*this, 0); ++i) {
                for (SizeType j = 0; j < extent(*this, 1); ++j) {
                    for (SizeType k = 0; k < extent(*this, 2); ++k) {
                        for (SizeType l = 0; l < extent(*this, 3); ++l) {
                            const Scalar val = 0.5 * ((i == k) && (j == l)) + 0.5 * ((i == l) && (j == k));
                            set(i, j, k, l, val);
                        }
                    }
                }
            }
        }

    private:
        ArrayView4D view_{};

        UTOPIA_FUNCTION TensorView(const TensorView &) { UTOPIA_DEVICE_ASSERT(false); }
    };

}  // namespace utopia

#endif  // UTOPIA_KOKKOS_TENSOR_VIEW_4_HPP
