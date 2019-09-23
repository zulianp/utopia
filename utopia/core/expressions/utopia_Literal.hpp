//
// Created by Patrick Zulian on 18/05/15.
//

#ifndef utopia_utopia_LITERAL_HPP
#define utopia_utopia_LITERAL_HPP

#include <iostream>
#include <cmath>
#include "utopia_Expression.hpp"
#include "utopia_BLAS_Operands.hpp"
#include "utopia_Normed.hpp"

namespace utopia {
    template<typename _Scalar>
    class Number : public Expression< Number<_Scalar> >, 
                   public BLAS1Tensor< Number<_Scalar> >,
                   public Normed<_Scalar> {
    public:
        static const int Order = 0;
        // static const int  StoreAs = UTOPIA_BY_VALUE;

        enum {
             StoreAs = UTOPIA_BY_VALUE
        };

        typedef _Scalar Scalar;
        // inline operator Scalar() const
        // {
        //     return value_;
        // }

        inline operator Scalar &()
        {
            return value_;
        }

        inline constexpr operator const Scalar &() const
        {
            return value_;
        }

        inline Number &operator+=(const Scalar &other)
        {
            value_ += other;
            return *this;
        }

        inline Number &operator-=(const Scalar &other)
        {
            value_ -= other;
            return *this;
        }

        inline Number &operator*=(const Scalar &other)
        {
            value_ *= other;
            return *this;
        }

        inline Number &operator/=(const Scalar &other)
        {
            value_ /= other;
            return *this;
        }

        // template<typename TOther>
        // inline constexpr Number &operator=(const Number<TOther> &other)
        // {
        //     value_ = other.value_;
        //     return *this;
        // }

        // inline constexpr Number &operator=(const Scalar &other)
        // {
        //     value_ = other;
        //     return *this;
        // }

        inline Scalar get() const
        {
            return value_;
        }

        inline void set(const Scalar &value)
        {
            value_ = value;
        }
        inline constexpr Number(const Scalar &value = 0.0)  : value_(value)
        {}

        inline std::string get_class() const override
        {
            return "Number";
        }

        template<typename OtherScalar>
        constexpr Number(const Number<OtherScalar> &other)
        : value_(other.value_)
        {}

        template<class Derived>
        Number(const Expression<Derived> &expr)
        : value_(scalar_cast<Scalar>(expr.derived()))
        {}

        template<class Derived>
        inline void construct(const Expression<Derived> &expr)
        {
            value_ = scalar_cast<Scalar>(expr.derived());
        }

        template<class OtherScalar>
        inline void construct(const Number<OtherScalar> &expr)
        {
            value_ = scalar_cast<Scalar>(expr.value_);
        }

        inline void construct(const Number<Scalar> &expr)
        {
            value_ = expr.value_;
        }

        ///<Scalar>SWAP - swap x and y
        inline void swap(Number &x) override
        {
            std::swap(x.value_, value_);
        }

        ///<Scalar>SCAL - x = a*x
        inline void scale(const Scalar &a) override
        {
            value_ *= a;
        }

        ///<Scalar>COPY - copy x into y (this)
        inline void copy(const Number &x) override
        {
            value_ = x.value_;
        }

        ///<Scalar>AXPY - y = a*x + y
        inline void axpy(const Scalar &a, const Number &x) override
        {
            value_ += a * x.value_;
        }

        ///<Scalar>DOT - dot product
        inline Scalar dot(const Number &other) const override
        {
            return value_ * other.value_;
        }

        ///<Scalar>NRM2 - Euclidean norm
        inline Scalar norm2() const override
        {
            return std::abs(value_);
        }

        ///<Scalar>ASUM - sum of absolute values
        inline Scalar norm1() const override
        {
            return std::abs(value_);
        }

        ///<Scalar>ASUM - sum of absolute values
        inline Scalar norm_infty() const override
        {
            return std::abs(value_);
        }

        ///I<Scalar>AMAX - index of max abs value
        inline unsigned short amax() const override
        {
            return 0;
        }


    private:
        Scalar value_;
    };

    template<typename T>
    Size size(const Number<T> &)
    {
        return { 1 };
    }

    template<typename Scalar_>
    class Traits< Number<Scalar_>> : public Traits<Scalar_> {
    public:
        using SizeType = unsigned short;
    };



   template<typename T>
   class IsScalar
   {
      public: enum { value = false };
   };


   template<typename T>
   class IsScalar< Number<T> >
   {
     public:  enum { value = true };
   };



    template<typename T>
    inline constexpr bool is_scalar()
    {
        return IsScalar<T>::value;
    }


    template<typename _Scalar, _Scalar Value = 0>
    class Literal : public Expression< Literal<_Scalar> > {
    public:
        typedef _Scalar Scalar;

        static const int Order = 0;

        inline constexpr operator Scalar() const
        {
            return Value;
        }

        inline std::string get_class() const
        {
            return "Literal(TODO)";
        }
    };

    template<typename T>
    void disp(const Number<T> &num, std::ostream &os = std::cout)
    {
        os << static_cast<T>(num) << std::endl;
    }
}
#endif //utopia_utopia_LITERAL_HPP
