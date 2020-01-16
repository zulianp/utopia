#ifndef UTOPIA_DEVICE_OUTER_PRODUCT_HPP
#define UTOPIA_DEVICE_OUTER_PRODUCT_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_ViewForwardDeclarations.hpp"
#include "utopia_StoreAs.hpp"

namespace utopia {

    template<class Left, class Right>
    class DeviceOuterProduct : public DeviceExpression<DeviceOuterProduct<Left, Right>> {
    public:
        using SizeType = typename Traits<Left>::SizeType;
        using Scalar   = typename Traits<Left>::Scalar;

        UTOPIA_INLINE_FUNCTION DeviceOuterProduct(const Left &left, const Right &right)
        : left_(left), right_(right)
        {}

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i, const SizeType &j) const
        {
            return left_[i] * right_[j];
        }

        UTOPIA_INLINE_FUNCTION SizeType rows() const
        {
            return left_.size();
        }

        UTOPIA_INLINE_FUNCTION SizeType cols() const
        {
            return right_.size();
        }

        friend UTOPIA_INLINE_FUNCTION SizeType rows(const DeviceOuterProduct &expr)
        {
            return expr.rows();
        }

        friend UTOPIA_INLINE_FUNCTION SizeType cols(const DeviceOuterProduct &expr)
        {
            return expr.cols();
        }

    private:
        UTOPIA_STORE_CONST(Left)  left_;
        UTOPIA_STORE_CONST(Right) right_;
    };

    template<class Left, class Right>
    class Traits<DeviceOuterProduct<Left, Right>> : public Traits<Left> {
    public:
        static const int Order = 2;
    };


}

#endif //UTOPIA_DEVICE_OUTER_PRODUCT_HPP
