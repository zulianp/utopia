#ifndef UTOPIA_VIEW_BINARY_HPP
#define UTOPIA_VIEW_BINARY_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_ViewForwardDeclarations.hpp"

#include "utopia_InlineEval.hpp"
#include "utopia_DeviceOp.hpp"
#include "utopia_DeviceExpression.hpp"

namespace utopia {

    template<class Left, class Right, class Op>
    class DeviceBinary : public DeviceExpression<DeviceBinary<Left, Right, Op>> {
    public:
        using SizeType = typename Traits<Left>::SizeType;
        using Scalar   = typename Traits<Left>::Scalar;

        UTOPIA_INLINE_FUNCTION DeviceBinary(const Left &left, const Right &right)
        : left_(left), right_(right)
        {
            // UTOPIA_DEVICE_ASSERT(left.size() == right.size());
        }

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i, const SizeType &j) const
        {
            return DeviceOp<Scalar, Op>::apply(left_(i, j), right_(i, j));
        }

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i) const
        {
            return DeviceOp<Scalar, Op>::apply(left_(i), right_(i));
        }

        inline std::string get_class() const override
        {
            return std::string("DeviceBinary<") + left_.get_class() + ", " + right_.get_class() + ", " + GetClass<Op>() + ">";
        }

        inline SizeType size() const
        {
            return left_.size();
        }

        UTOPIA_INLINE_FUNCTION const Left &left() const
        {
            return left_;
        }

        UTOPIA_INLINE_FUNCTION const Right &right() const
        {
            return right_;
        }

        UTOPIA_INLINE_FUNCTION SizeType rows() const
        {
            return right_.rows();
        }

        UTOPIA_INLINE_FUNCTION SizeType cols() const
        {
            return right_.cols();
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

        UTOPIA_INLINE_FUNCTION const Left &left() const
        {
            return left_;
        }

        UTOPIA_INLINE_FUNCTION const Right &right() const
        {
            return right_;
        }

        UTOPIA_INLINE_FUNCTION SizeType rows() const
        {
            return right_.rows();
        }

        UTOPIA_INLINE_FUNCTION SizeType cols() const
        {
            return right_.cols();
        }

    private:
        const Left left_;
        UTOPIA_STORE_CONST(Right) right_;
    };

    template<class Left, class Right, class Operation>
    class Traits< DeviceBinary<Left, Right, Operation> > : public Traits< typename MostDescriptive<Left, Right>::Type > {};

    template<class Left, class Right, class Operation>
    UTOPIA_INLINE_FUNCTION
    typename Traits<DeviceBinary<Left, Right, Operation>>::SizeType
    rows(const DeviceBinary<Left, Right, Operation> &expr)
    {
        return expr.rows();
    }

    template<class Left, class Right, class Operation>
    UTOPIA_INLINE_FUNCTION
    typename Traits<DeviceBinary<Left, Right, Operation> >::SizeType
    cols(const DeviceBinary<Left, Right, Operation> &expr)
    {
        return expr.cols();
    }
}

#endif
