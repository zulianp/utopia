#ifndef UTOPIA_DEVICE_TENSOR_PRODUCT_HPP
#define UTOPIA_DEVICE_TENSOR_PRODUCT_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_ViewForwardDeclarations.hpp"
#include "utopia_StoreAs.hpp"

namespace utopia {

    template<class Left, class Right, int... >
    class DeviceTensorProduct {};

    template<typename SizeType, int Index, int NIndices>
    class TensorIndexSelector {};


    template<typename SizeType>
    class TensorIndexSelector<SizeType, 0, 4> {
    public:
        static constexpr UTOPIA_INLINE_FUNCTION SizeType select(
            const SizeType &i,
            const SizeType &,
            const SizeType &,
            const SizeType &
        )
        {
            return i;
        }
    };

    template<typename SizeType>
    class TensorIndexSelector<SizeType, 1, 4> {
    public:
        static constexpr UTOPIA_INLINE_FUNCTION SizeType select(
            const SizeType &,
            const SizeType &j,
            const SizeType &,
            const SizeType &
        )
        {
            return j;
        }
    };

    template<typename SizeType>
    class TensorIndexSelector<SizeType, 2, 4> {
    public:
        static constexpr UTOPIA_INLINE_FUNCTION SizeType select(
            const SizeType &,
            const SizeType &,
            const SizeType &k,
            const SizeType &
        )
        {
            return k;
        }
    };

    template<typename SizeType>
    class TensorIndexSelector<SizeType, 3, 4> {
    public:
        static constexpr UTOPIA_INLINE_FUNCTION SizeType select(
            const SizeType &,
            const SizeType &,
            const SizeType &,
            const SizeType &l
        )
        {
            return l;
        }
    };

    template<class Left, class Right, int I, int J, int K, int L>
    class DeviceTensorProduct<Left, Right, I, J, K, L> : public DeviceExpression<DeviceTensorProduct<Left, Right, I, J, K, L>> {
    public:
        using SizeType = typename Traits<Left>::SizeType;
        using Scalar   = typename Traits<Left>::Scalar;

        using ISelector = utopia::TensorIndexSelector<SizeType, I, 4>;
        using JSelector = utopia::TensorIndexSelector<SizeType, J, 4>;
        using KSelector = utopia::TensorIndexSelector<SizeType, K, 4>;
        using LSelector = utopia::TensorIndexSelector<SizeType, L, 4>;

        UTOPIA_INLINE_FUNCTION DeviceTensorProduct(const Left &left, const Right &right)
        : left_(left), right_(right)
        {}

        UTOPIA_INLINE_FUNCTION Scalar operator()(
            const SizeType &i,
            const SizeType &j,
            const SizeType &k,
            const SizeType &l) const
        {
            return left_(
                    ISelector::select(i, j, k, l),
                    JSelector::select(i, j, k, l)
                    ) *
                    right_(
                        KSelector::select(i, j, k, l),
                        LSelector::select(i, j, k, l)
                    );
        }

        // friend UTOPIA_INLINE_FUNCTION SizeType extent(const DeviceTensorProduct &that, const SizeType &i)
        // {

        // }


    private:
        UTOPIA_STORE_CONST(Left)  left_;
        UTOPIA_STORE_CONST(Right) right_;
    };


    template<class Left, class Right, int I, int J, int K, int L>
    class Traits<DeviceTensorProduct<Left, Right, I, J, K, L>> : public Traits<Right> {
    public:
        static const int Order = 4;
    };


    template<int...Ind, class Left, class Right>
    UTOPIA_INLINE_FUNCTION DeviceTensorProduct<Left, Right, Ind...> tensor_product(
        const DeviceExpression<Left> &left,
        const DeviceExpression<Right> &right)
    {
        return DeviceTensorProduct<Left, Right, Ind...>(left, right);
    }

}

#endif //UTOPIA_DEVICE_TENSOR_PRODUCT_HPP
