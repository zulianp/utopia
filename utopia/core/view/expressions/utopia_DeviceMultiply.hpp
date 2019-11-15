#ifndef UTOPIA_DEVICE_MULTIPLY_HPP
#define UTOPIA_DEVICE_MULTIPLY_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_ViewForwardDeclarations.hpp"
#include "utopia_StoreAs.hpp"
// #include "utopia_InlineEval.hpp"
// #include "utopia_DeviceOp.hpp"
// #include "utopia_DeviceExpression.hpp"

namespace utopia {

    template<class Left, class Right>
    class DeviceMultiply : public DeviceExpression<DeviceMultiply<Left, Right>> {
    public:
        using SizeType = typename Traits<Left>::SizeType;
        using Scalar   = typename Traits<Left>::Scalar;

        UTOPIA_INLINE_FUNCTION DeviceMultiply(const Left &left, const Right &right) 
        : left_(left), right_(right)
        {
            // UTOPIA_DEVICE_ASSERT(left.cols() == right.size());
        }

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i, const SizeType &j) const
        {
            const SizeType cols = left_.cols();

            Scalar ret = 0.0;
            for(SizeType k = 0; k < cols; ++k) {
                ret += left_(i, k) * right_(k, j);
            }

            return ret;
        }

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i) const
        {
            const SizeType cols = left_.cols();

            Scalar ret = 0.0;
            for(SizeType j = 0; j < cols; ++j) {
                ret += left_(i, j) * right_(j);
            }

            return ret;
        }

        inline std::string get_class() const override
        {
            return "DeviceMultiply<" + left_.get_class() + ", " + right_.get_class() + ">";
        }

    private:
        UTOPIA_STORE_CONST(Left)  left_;
        UTOPIA_STORE_CONST(Right) right_;

    };

    template<class Left, class Right>
    class Traits< DeviceMultiply<Left, Right> > : public Traits<typename ChooseType<Left, Right, Right>::Type > {};
}

#endif
