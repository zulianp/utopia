#ifndef UTOPIA_DEVICE_IN_PLACE_HPP
#define UTOPIA_DEVICE_IN_PLACE_HPP

namespace utopia {
    template<class Left, class Right, class Op, int Order = Traits<Left>::Order>
    class DeviceInPlace {};

    template<class Left, class Right, class Op>
    class DeviceInPlace<Left, Right, Op, 1> {
    public:
        using Scalar = typename Traits<Left>::Scalar;
        using SizeType = typename Traits<Left>::SizeType;

        UTOPIA_INLINE_FUNCTION static void apply(Left &left, const Right &right)
        {
            const SizeType n = left.size();
            for(SizeType i = 0; i < n; ++i) {
                left(i) = DeviceOp<Scalar, Op>::apply(left(i), right(i));
            }
        }
    };

    template<class Left, class Right, class Op>
    class DeviceInPlace<Left, Right, Op, 2> {
    public:
        using Scalar   = typename Traits<Left>::Scalar;
        using SizeType = typename Traits<Left>::SizeType;

        UTOPIA_INLINE_FUNCTION static void apply(Left &left, const Right &right)
        {
            const SizeType rows = left.rows();
            const SizeType cols = left.cols();

            for(SizeType i = 0; i < rows; ++i) {
                for(SizeType j = 0; j < cols; ++j) {
                    left(i, j) = DeviceOp<Scalar, Op>::apply(left(i, j), right(i, j));
                }
            }
        }
    };

    template<class Left, class Right, class Op>
    class DeviceInPlace<Left, Right, Op, 0> {
    public:
        using Scalar = typename Traits<Left>::Scalar;

        UTOPIA_INLINE_FUNCTION void apply(Left &left, const Right &right)
        {
            left = DeviceOp<Scalar, Op>::apply(left, right);
        }
    };
}

#endif //UTOPIA_DEVICE_IN_PLACE_HPP
