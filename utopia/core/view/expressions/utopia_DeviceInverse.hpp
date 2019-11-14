#ifndef UTOPIA_DEVICE_INVERSE_HPP
#define UTOPIA_DEVICE_INVERSE_HPP


namespace utopia {

    template<class Left, class Right>
    class DeviceInverse {
    public:
        using Scalar   = typename Traits<Left>::Scalar;
        using SizeType = typename Traits<Left>::SizeType;

        UTOPIA_INLINE_FUNCTION static bool apply(Left &left, const Right &right)
        {

            return false;
        }
    };

}

#endif //UTOPIA_DEVICE_INVERSE_HPP
