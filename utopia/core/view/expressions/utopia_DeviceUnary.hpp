#ifndef UTOPIA_DEVICE_UNARY_HPP
#define UTOPIA_DEVICE_UNARY_HPP

#include "utopia_Base.hpp"
#include "utopia_DeviceExpression.hpp"
#include "utopia_DeviceOp.hpp"
#include "utopia_InlineEval.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template <class InnerExpr, class Op>
    class DeviceUnary : public DeviceExpression<DeviceUnary<InnerExpr, Op>> {
    public:
        using SizeType = typename Traits<InnerExpr>::SizeType;
        using Scalar = typename Traits<InnerExpr>::Scalar;

        UTOPIA_INLINE_FUNCTION DeviceUnary(const InnerExpr &expr) : expr_(expr) {}

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i, const SizeType &j) const {
            return DeviceOp<Scalar, Op>::apply(expr_(i, j));
        }

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i) const {
            return DeviceOp<Scalar, Op>::apply(expr_(i));
        }

        UTOPIA_INLINE_FUNCTION SizeType rows() const { return expr_.rows(); }

        UTOPIA_INLINE_FUNCTION SizeType cols() const { return expr_.cols(); }

        UTOPIA_INLINE_FUNCTION SizeType size() const { return expr_.size(); }

    private:
        UTOPIA_STORE_CONST(InnerExpr) expr_;
    };

    template <class InnerExpr, class Op>
    class Traits<DeviceUnary<InnerExpr, Op>> : public Traits<InnerExpr> {};
}  // namespace utopia

#endif
