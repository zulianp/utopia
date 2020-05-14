#ifndef UTOPIA_DEVICE_NUMBER_HPP
#define UTOPIA_DEVICE_NUMBER_HPP

#include "utopia_ViewForwardDeclarations.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Literal.hpp"

namespace utopia {

    template<typename _Scalar>
    class DeviceNumber : public DeviceExpression< DeviceNumber<_Scalar> > {
    public:
        static const int Order = 0;
        // static const int StoreAs = UTOPIA_BY_VALUE;

        enum {
             StoreAs = UTOPIA_BY_VALUE
        };

        typedef _Scalar Scalar;

        UTOPIA_INLINE_FUNCTION operator Scalar &()
        {
            return value_;
        }

        UTOPIA_INLINE_FUNCTION constexpr operator const Scalar &() const
        {
            return value_;
        }

        UTOPIA_INLINE_FUNCTION DeviceNumber &operator+=(const Scalar &other)
        {
            value_ += other;
            return *this;
        }

        UTOPIA_INLINE_FUNCTION DeviceNumber &operator-=(const Scalar &other)
        {
            value_ -= other;
            return *this;
        }

        UTOPIA_INLINE_FUNCTION DeviceNumber &operator*=(const Scalar &other)
        {
            value_ *= other;
            return *this;
        }

        UTOPIA_INLINE_FUNCTION DeviceNumber &operator/=(const Scalar &other)
        {
            value_ /= other;
            return *this;
        }

        UTOPIA_INLINE_FUNCTION bool is_alias(const DeviceNumber &other) const
        {
            return this == &other;
        }

        UTOPIA_INLINE_FUNCTION Scalar get() const
        {
            return value_;
        }

        UTOPIA_INLINE_FUNCTION void set(const Scalar &value)
        {
            value_ = value;
        }
        UTOPIA_INLINE_FUNCTION constexpr DeviceNumber(const Scalar &value = Zero<Scalar>::value())  : value_(value)
        {}

        inline std::string get_class() const override /*override*/
        {
            return "DeviceNumber";
        }

        template<typename OtherScalar>
        UTOPIA_INLINE_FUNCTION constexpr DeviceNumber(const DeviceNumber<OtherScalar> &other)
        : value_(other.value_)
        {}

        template<typename OtherScalar>
        UTOPIA_INLINE_FUNCTION constexpr DeviceNumber(DeviceNumber<OtherScalar> &&other)
        : value_(std::move(other.value_))
        {}

        // template<class Derived>
        // DeviceNumber(const DeviceExpression<Derived> &expr)
        // : value_(scalar_cast<Scalar>(expr.derived()))
        // {}

        template<class Derived>
        UTOPIA_INLINE_FUNCTION void construct(const DeviceExpression<Derived> &expr)
        {
            value_ = scalar_cast<Scalar>(expr.derived());
        }

        template<class OtherScalar>
        UTOPIA_INLINE_FUNCTION void construct(const DeviceNumber<OtherScalar> &expr)
        {
            value_ = scalar_cast<Scalar>(expr.value_);
        }

        UTOPIA_INLINE_FUNCTION void construct(const DeviceNumber<Scalar> &expr)
        {
            value_ = expr.value_;
        }

        UTOPIA_INLINE_FUNCTION void assign(const DeviceNumber<Scalar> &expr)
        {
            value_ = expr.value_;
        }

        ///<Scalar>SWAP - swap x and y
        UTOPIA_INLINE_FUNCTION void swap(DeviceNumber &x) /*override*/
        {
            std::swap(x.value_, value_);
        }

        ///<Scalar>SCAL - x = a*x
        UTOPIA_INLINE_FUNCTION void scale(const Scalar &a) /*override*/
        {
            value_ *= a;
        }

        ///<Scalar>COPY - copy x into y (this)
        UTOPIA_INLINE_FUNCTION void copy(const DeviceNumber &x) /*override*/
        {
            value_ = x.value_;
        }

        ///<Scalar>AXPY - y = a*x + y
        UTOPIA_INLINE_FUNCTION void axpy(const Scalar &a, const DeviceNumber &x) /*override*/
        {
            value_ += a * x.value_;
        }

        ///<Scalar>DOT - dot product
        UTOPIA_INLINE_FUNCTION Scalar dot(const DeviceNumber &other) const /*override*/
        {
            return value_ * other.value_;
        }

        ///<Scalar>NRM2 - Euclidean norm
        UTOPIA_INLINE_FUNCTION Scalar norm2() const /*override*/
        {
            return DeviceOp<Scalar, Abs>::apply(value_);
        }

        ///<Scalar>ASUM - sum of absolute values
        UTOPIA_INLINE_FUNCTION Scalar norm1() const /*override*/
        {
            return DeviceOp<Scalar, Abs>::apply(value_);
        }

        UTOPIA_INLINE_FUNCTION Scalar norm_infty() const /*override*/
        {
            return DeviceOp<Scalar, Abs>::apply(value_);
        }

        ///I<Scalar>AMAX - index of max abs value
        UTOPIA_INLINE_FUNCTION unsigned short amax() const ///*override*/
        {
            return 0;
        }

        UTOPIA_INLINE_FUNCTION static constexpr unsigned short size() { return 1; }
        UTOPIA_INLINE_FUNCTION static constexpr unsigned short rows() { return 1; }
        UTOPIA_INLINE_FUNCTION static constexpr unsigned short cols() { return 1; }

    private:
        Scalar value_;
    };

    template<typename Scalar>
    class Traits< DeviceNumber<Scalar>> : public Traits<Scalar> {
    public:
        static const int Order = 0;
    };

    template<class Left, class Right>
    class MostDescriptive<Left, DeviceNumber<Right> > {
    public:
        using Type = typename Unwrap<Left>::Type;
    };

    template<class Left, class Right>
    class MostDescriptive<DeviceNumber<Left>, Right > {
    public:
        using Type = typename Unwrap<Right>::Type;
    };

    template<typename Left, typename Right>
    class MostDescriptive<DeviceNumber<Left>, DeviceNumber<Right> > {
    public:
        // typedef decltype(Left(0) + Right(0)) Type;
        using Type = utopia::DeviceNumber<decltype(Left(0) + Right(0))>;
    };

    template<typename Scalar>
    UTOPIA_INLINE_FUNCTION constexpr int size(const DeviceNumber<Scalar> &)
    {
       return 1;
    }

    template<typename Scalar>
    UTOPIA_INLINE_FUNCTION constexpr int rows(const DeviceNumber<Scalar> &)
    {
       return 1;
    }

    template<typename Scalar>
    UTOPIA_INLINE_FUNCTION constexpr int cols(const DeviceNumber<Scalar> &)
    {
       return 1;
    }

}

#endif //UTOPIA_DEVICE_NUMBER_HPP
