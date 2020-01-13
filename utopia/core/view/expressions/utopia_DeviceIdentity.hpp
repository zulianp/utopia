#ifndef UTOPIA_DEVICEIDENTITY_HPP
#define UTOPIA_DEVICEIDENTITY_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_DeviceExpression.hpp"
#include "utopia_DeviceBinary.hpp"

namespace utopia {

    template<typename Scalar_>
    class DeviceIdentity : public DeviceExpression<DeviceIdentity<Scalar_>> {
    public:
        using Scalar = Scalar_;

        template<typename SizeType>
        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i, const SizeType &j) const
        {
            return static_cast<Scalar>((i == j) ? 1 : 0);
        }

        // UTOPIA_INLINE_FUNCTION operator Scalar() const
        // {
        //     return 1.0;
        // }
    };

    template<typename T>
    class HasSize<DeviceIdentity<T>> {
    public:
        static const int value = 0;
    };

    template<typename Scalar_>
    class Traits<DeviceIdentity<Scalar_>> {
    public:
        static const int Order = 2;
        using Scalar = Scalar_;
    };

    template<class Left, typename Scalar>
    class MostDescriptive<Left, DeviceIdentity<Scalar>> {
    public:
        using Type = typename Unwrap<Left>::Type;
    };

    template<typename Scalar, class Right>
    class MostDescriptive<DeviceIdentity<Scalar>, Right> {
    public:
        using Type = typename Unwrap<Right>::Type;
    };

    template<typename Left, typename Right>
    class MostDescriptive<DeviceIdentity<Left>, DeviceIdentity<Right>> {
    public:
        using Type = DeviceIdentity<decltype(Left(0.0) * Right(0.0))>;
    };


    template<typename Left, typename Right>
    class MostDescriptive<DeviceNumber<Left>, DeviceIdentity<Right>> {
    public:
        using Type = DeviceIdentity<decltype(Left(0.0) * Right(0.0))>;
    };


    namespace device {

        template<typename Scalar>
        inline constexpr DeviceIdentity<Scalar> identity()
        {
            return DeviceIdentity<Scalar>();
        }
    }


    template<class Left, class Right, class Op>
    class DeviceBinary<DeviceNumber<Left>, DeviceIdentity<Right>, Op> :
        public DeviceExpression<DeviceBinary<DeviceNumber<Left>, DeviceIdentity<Right>, Op> >{
    public:
        using Scalar   = typename Traits<Right>::Scalar;


        UTOPIA_INLINE_FUNCTION DeviceBinary(const DeviceNumber<Left> &left, const DeviceIdentity<Right> &right) : left_(left), right_(right) {}

        template<typename SizeType>
        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i, const SizeType &j) const
        {
            return DeviceOp<Scalar, Op>::apply(left_, right_(i, j));
        }

        template<typename SizeType>
        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i) const
        {
            return DeviceOp<Scalar, Op>::apply(left_, right_(i));
        }

        UTOPIA_INLINE_FUNCTION const Left &left() const
        {
            return left_;
        }

        UTOPIA_INLINE_FUNCTION const DeviceIdentity<Right> &right() const
        {
            return right_;
        }

    private:
        Left left_;
        DeviceIdentity<Right> right_;
    };


    template<class Left, class Right, class Op>
    class HasSize< DeviceBinary<DeviceNumber<Left>, DeviceIdentity<Right>, Op> > {
    public:
        static const int value = 0;
    };

}

#endif //UTOPIA_DEVICEIDENTITY_HPP