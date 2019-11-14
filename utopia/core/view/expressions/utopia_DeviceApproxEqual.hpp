#ifndef UTOPIA_DEVICE_APPROX_EQUAL_HPP
#define UTOPIA_DEVICE_APPROX_EQUAL_HPP

namespace utopia {
    template<class Left, class Right, int Order = Traits<Left>::Order>
    class DeviceApproxEqual {};

    template<class Left, class Right>
    class DeviceApproxEqual<Left, Right, 1> {
    public:
        using Scalar = typename Traits<Left>::Scalar;
        using SizeType = typename Traits<Left>::SizeType;

        UTOPIA_INLINE_FUNCTION static bool apply(const Left &left, const Right &right, const Scalar &tol)
        {
            const SizeType n = left.size();
            bool ret = n != 0;
            for(SizeType i = 0; i < n; ++i) {
                if(!device::approxeq(left(i), right(i), tol)) {
                    ret = false;
                    break;
                }
            }

            return ret;
        }
    };

    template<class Left, class Right>
    class DeviceApproxEqual<Left, Right, 2> {
    public:
        using Scalar   = typename Traits<Left>::Scalar;
        using SizeType = typename Traits<Left>::SizeType;

        UTOPIA_INLINE_FUNCTION static bool apply(const Left &left, const Right &right, const Scalar &tol)
        {
            const SizeType rows = left.rows();
            const SizeType cols = left.cols();

            bool ret = rows != 0;
            
            for(SizeType i = 0; i < rows; ++i) {
                for(SizeType j = 0; j < cols; ++j) {
                    if(!device::approxeq(left(i, j), right(i, j), tol)) {
                        ret = false;
                        break;
                    }
                }
            }

            return ret;
        }
    };

    template<class Left, class Right>
    class DeviceApproxEqual<Left, Right, 0> {
    public:
        using Scalar = typename Traits<Left>::Scalar;

        UTOPIA_INLINE_FUNCTION bool apply(const Left &left, const Right &right, const Scalar &tol)
        {
            return device::approxeq(left, right, tol);
        }
    };

    template<class Left, class Right>
    class Traits< DeviceApproxEqual<Left, Right> > : public Traits< typename MostDescriptive<Left, Right>::Type > {};
}

#endif
