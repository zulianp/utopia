#ifndef UTOPIA_DEVICE_DOT_HPP
#define UTOPIA_DEVICE_DOT_HPP

namespace utopia {

    template<class Left, class Right, int Order = Traits<Left>::Order>
    class DeviceDot {};

    template<class Left, class Right>
    class DeviceDot<Left, Right, 1> {
    public:
        using Scalar   = typename Traits<Right>::Scalar;
        using SizeType = typename Traits<Right>::SizeType;

        UTOPIA_INLINE_FUNCTION static Scalar apply(
            const Left &left,
            const Right &right)
        {
            Scalar ret = 0.0;

            const SizeType n = left.size();
            for(SizeType i = 0; i < n; ++i) {
                ret += left(i) * right(i);
            }

            return ret;
        }
    };


    template<class Left, class Right>
    class DeviceDot<Left, Right, 2> {
    public:
        using Scalar   = typename Traits<Right>::Scalar;
        using SizeType = typename Traits<Right>::SizeType;

        UTOPIA_INLINE_FUNCTION static Scalar apply(
            const Left &left,
            const Right &right)
        {
            Scalar ret = 0.0;

            const SizeType rows = left.rows();
            const SizeType cols = left.cols();

            for(SizeType i = 0; i < rows; ++i) {
                for(SizeType j = 0; j < cols; ++j) {
                    ret += left(i, j) * right(i, j);
                }
            }

            return ret;
        }
    };
  
}

#endif //UTOPIA_DEVICE_DOT_HPP
