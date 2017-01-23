//
// Created by Patrick Zulian on 18/05/15.
//

#ifndef utopia_utopia_LITERAL_HPP
#define utopia_utopia_LITERAL_HPP

#include <iostream>
#include "utopia_Expression.hpp"

namespace utopia {
    template<typename _Scalar>
    class Number : public Expression< Number<_Scalar> > {
    public:
        enum {
         Order = 0
        };

         enum {
            StoreAs = UTOPIA_BY_VALUE
        };


        typedef _Scalar Scalar;
        // inline operator Scalar() const
        // {
        //     return _value;
        // }

        inline operator Scalar &()
        {
            return _value;
        }

        inline operator const Scalar &() const
        {
            return _value;
        }


        inline Number(const Scalar &value)  : _value(value)
        {}

        std::string getClass() const
        {
            return "Number";
        }

        template<typename OtherScalar>
        Number(const Number<OtherScalar> &other)
        : _value(other._value)
        {}

        template<class Derived>
        Number(const Expression<Derived> &expr)
        : _value(scalar_cast<Scalar>(expr))
        {}

    private:
        Scalar _value;
    };

    template<typename T>
    Size size(const Number<T> &)
    {
        return { 1 };
    }



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

        enum {
            Order = 0
        };

        inline constexpr operator Scalar() const
        {
            return Value;
        }

        inline std::string getClass() const
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
