#ifndef UTOPIA_DEVICE_TRANSPOSE_HPP
#define UTOPIA_DEVICE_TRANSPOSE_HPP

#include "utopia_Traits.hpp"
#include "utopia_ViewForwardDeclarations.hpp"

#include "utopia_MatrixView.hpp"

namespace utopia {

    template <class Expr>
    class DeviceTranspose : public DeviceExpression<DeviceTranspose<Expr>> {
    public:
        using Scalar = typename Traits<Expr>::Scalar;
        using SizeType = typename Traits<Expr>::SizeType;

        UTOPIA_FUNCTION DeviceTranspose(const Expr &expr) : expr_(expr) {}

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i, const SizeType &j) const {
            return DeviceOp<Scalar, Conjugate>::apply(expr_(j, i));
        }

        UTOPIA_INLINE_FUNCTION const Expr &expr() const { return expr_; }

        UTOPIA_INLINE_FUNCTION SizeType rows() const { return expr_.cols(); }

        UTOPIA_INLINE_FUNCTION SizeType cols() const { return expr_.rows(); }

    private:
        UTOPIA_STORE_CONST(Expr) expr_;
    };

    template <class InnerExpr>
    class Traits<DeviceTranspose<InnerExpr>> : public Traits<InnerExpr> {};

    template <class Expr>
    UTOPIA_INLINE_FUNCTION typename Traits<Expr>::SizeType rows(const DeviceTranspose<Expr> &t) {
        return utopia::cols(t.expr());
    }

    template <class Expr>
    UTOPIA_INLINE_FUNCTION typename Traits<Expr>::SizeType cols(const DeviceTranspose<Expr> &t) {
        return utopia::rows(t.expr());
    }

    template <class Derived>
    UTOPIA_INLINE_FUNCTION DeviceTranspose<Derived> transpose(const DeviceExpression<Derived> &expr) {
        return DeviceTranspose<Derived>(expr.derived());
    }

    template <class Derived>
    UTOPIA_INLINE_FUNCTION UTOPIA_STORE_CONST(Derived) transpose(const DeviceTranspose<Derived> &expr) {
        return expr.expr();
    }

}  // namespace utopia

#endif  // UTOPIA_DEVICE_TRANSPOSE_HPP
