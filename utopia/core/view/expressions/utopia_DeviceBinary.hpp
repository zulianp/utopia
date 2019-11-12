#ifndef UTOPIA_VIEW_BINARY_HPP
#define UTOPIA_VIEW_BINARY_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_InlineEval.hpp"
#include "utopia_DeviceOp.hpp"
#include "utopia_DeviceExpression.hpp"

namespace utopia {

    template<class Left, class Right, class Op>
    class DeviceBinary :
        public DeviceExpression<DeviceBinary<Left, Right, Op>>{
    public:
        using SizeType = typename Traits<Left>::SizeType;
        using Scalar   = typename Traits<Left>::Scalar;

        UTOPIA_INLINE_FUNCTION DeviceBinary(const Left &left, const Right &right) : left_(left), right_(right) {}

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i, const SizeType &j) const
        {
            return DeviceOp<Scalar, Op>::apply(left_(i, j), right_(i, j));
        }

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i) const
        {
            return DeviceOp<Scalar, Op>::apply(left_(i), right_(i));
        }

    private:
        UTOPIA_STORE_CONST(Left)  left_;
        UTOPIA_STORE_CONST(Right) right_;

    };

    template<class Left, class Right, class Op>
    class DeviceBinary<DeviceNumber<Left>, Right, Op> :
        public DeviceExpression<DeviceBinary<DeviceNumber<Left>, Right, Op> >{
    public:
        using SizeType = typename Traits<Right>::SizeType;
        using Scalar   = typename Traits<Right>::Scalar;

        UTOPIA_INLINE_FUNCTION DeviceBinary(const DeviceNumber<Left> &left, const Right &right) : left_(left), right_(right) {}

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i, const SizeType &j) const
        {
            return DeviceOp<Scalar, Op>::apply(left_, right_(i, j));
        }

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i) const
        {
            return DeviceOp<Scalar, Op>::apply(left_, right_(i));
        }

    private:
        const Left left_;
        UTOPIA_STORE_CONST(Right) right_;
    };
}

#endif