#ifndef UTOPIA_DEVICE_TENSOR_CONTRACTION_HPP
#define UTOPIA_DEVICE_TENSOR_CONTRACTION_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_ViewForwardDeclarations.hpp"
#include "utopia_StoreAs.hpp"

namespace utopia {

    template<class Left, class Right, int LeftOrder = Traits<Left>::Order, int RightOrder = Traits<Right>::Order>
    class DeviceTensorContraction {};

    template<class Left, class Right>
    class DeviceTensorContraction<Left, Right, 4, 2> : public DeviceExpression<DeviceTensorContraction<Left, Right, 4, 2>> {
    public:
        using SizeType = typename Traits<Left>::SizeType;
        using Scalar   = typename Traits<Left>::Scalar;

        UTOPIA_INLINE_FUNCTION DeviceTensorContraction(const Left &left, const Right &right)
        : left_(left), right_(right)
        {}

        UTOPIA_INLINE_FUNCTION Scalar operator()(
            const SizeType &i,
            const SizeType &j) const
        {
            const auto rows = right_.rows();
            const auto cols = right_.cols();

            Scalar ret = 0.0;

            for(SizeType k = 0; k < rows; ++k) {
                for(SizeType l = 0; l < cols; ++l) {
                    ret += left_(i, j, k, l) * right_(k, l);
                }
            }

            return ret;
        }

    private:
        UTOPIA_STORE_CONST(Left)  left_;
        UTOPIA_STORE_CONST(Right) right_;

    };


    template<class Left, class Right>
    class Traits< DeviceTensorContraction<Left, Right, 4, 2> > : public Traits<Right> {};


    template<class Left, class Right>
    UTOPIA_INLINE_FUNCTION DeviceTensorContraction<Left, Right> contraction(
        const DeviceExpression<Left> &left,
        const DeviceExpression<Right> &right)
    {
        return DeviceTensorContraction<Left, Right>(left, right);
    }
}

#endif //UTOPIA_TENSOR_CONTRACTION_HPP
