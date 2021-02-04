#ifndef UTOPIA_DEVICE_CROSS_PRODUCT_HPP
#define UTOPIA_DEVICE_CROSS_PRODUCT_HPP

#include "utopia_Base.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"
#include "utopia_ViewForwardDeclarations.hpp"

#include "utopia_DeviceExpression.hpp"
#include "utopia_DeviceOp.hpp"
#include "utopia_InlineEval.hpp"

namespace utopia {

    template <class Left, class Right>
    class DeviceCrossProduct : public DeviceExpression<DeviceCrossProduct<Left, Right>> {
    public:
        using SizeType = typename Traits<Left>::SizeType;
        using Scalar = typename Traits<Left>::Scalar;

        UTOPIA_INLINE_FUNCTION DeviceCrossProduct(const Left &left, const Right &right) : left_(left), right_(right) {
            UTOPIA_DEVICE_ASSERT(left.size() == right.size());
            UTOPIA_DEVICE_ASSERT(left.size() == 3);
        }

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i) const {
            switch (i) {
                case 0: {
                    return (left_(1) * right_(2)) - (left_(2) * right_(1));
                }

                case 1: {
                    return (left_(2) * right_(0)) - (left_(0) * right_(2));
                }

                case 2: {
                    return (left_(0) * right_(1)) - (left_(1) * right_(0));
                }

                default: {
                    return 0;
                }
            }
        }

        inline std::string get_class() const override {
            return std::string("DeviceCrossProduct<") + left_.get_class() + ", " + right_.get_class() + ">";
        }

        UTOPIA_INLINE_FUNCTION const Left &left() const { return left_; }

        UTOPIA_INLINE_FUNCTION const Right &right() const { return right_; }

        UTOPIA_INLINE_FUNCTION static constexpr SizeType size() { return 3; }

    private:
        UTOPIA_STORE_CONST(Left) left_;
        UTOPIA_STORE_CONST(Right) right_;
    };

    template <class Left, class Right>
    class Traits<DeviceCrossProduct<Left, Right>> : public Traits<typename MostDescriptive<Left, Right>::Type> {};

}  // namespace utopia

#endif  // UTOPIA_DEVICE_CROSS_PRODUCT_HPP