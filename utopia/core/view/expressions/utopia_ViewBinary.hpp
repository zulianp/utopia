#ifndef UTOPIA_VIEW_BINARY_HPP
#define UTOPIA_VIEW_BINARY_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_InlineEval.hpp"
#include "utopia_DeviceOp.hpp"
#include "utopia_ViewExpression.hpp"

namespace utopia {

    template<class Left, class Right, class Op>
    class ViewBinary : public ViewExpression<ViewBinary<Left, Right, Op>>{
    public:
        using SizeType = typename Traits<Left>::SizeType;
        using Scalar = typename Traits<Left>::Scalar;

        UTOPIA_INLINE_FUNCTION ViewBinary(Left &left, const Right &right) : left_(left), right_(right) {}

        UTOPIA_INLINE_FUNCTION void apply(const SizeType &i, const SizeType &j) const
        {
            return DeviceOp<Scalar, Op>::apply(left_(i, j), right_(i, j));
        }

        UTOPIA_INLINE_FUNCTION void apply(const SizeType &i) const
        {
            return DeviceOp<Scalar, Op>::apply(left_(i), right_(i));
        }

    private:
        Left &left_;
        UTOPIA_STORE_CONST(Right) right_;

    };
}

#endif