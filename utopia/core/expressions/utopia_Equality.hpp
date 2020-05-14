#ifndef UTOPIA_EQUALITY_HPP
#define UTOPIA_EQUALITY_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Expression.hpp"


namespace utopia {
    template<class Left, class Right>
    class Equality : public Expression< Equality<Left, Right> > {
    public:

        static const int Order = Left::Order;

        using Scalar = typename Left::Scalar;

        Equality(const Left &left, const Right &right)
        : left_(left), right_(right)
        {}

        inline const Left &left() const
        {
            return left_;
        }

        inline const Right &right() const
        {
            return right_;
        }


        // inline Left &left()
        // {
        // 	return left_;
        // }

        // inline Right &right()
        // {
        // 	return right_;
        // }

        std::string get_class() const
        {
            return "Equality<" + left().get_class() + ", " + right().get_class() + ">";
        }

    private:
        UTOPIA_STORE_CONST(Left)  left_;
        UTOPIA_STORE_CONST(Right) right_;
    };

    template<class Left, class Right>
    class Traits< Equality<Left, Right> > : public Traits<Left> {};

    template<class Left, class Right>
    Size size(const Equality<Left, Right>  &expr)
    {
        return size(expr.left());
    }

    template<class Left, class Right>
    inline Equality<Left, Right> operator==(const Expression<Left> &left,
                                            const Expression<Right> &right)
    {
        return Equality<Left, Right>(left.derived(), right.derived());
    }

    template<class Left, class Right>
    inline Equality<Left, Number<Right> > operator==(const Expression<Left> &left,
                                                     const Number<Right> &right)
    {
        return Equality<Left, Number<Right> >(left.derived(), right);
    }

    ///If there is a scalar argument, then it is always place on the right
    template<class Left, class Right>
    inline Equality<Left, Number<Right> > operator==(const Number<Right> &right,
                                                     const Expression<Left> &left)
    {
        return Equality<Left, Number<Right> >(left.derived(), right);
    }
}

#endif //UTOPIA_EQUALITY_HPP

