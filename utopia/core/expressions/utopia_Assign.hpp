//
// Created by Patrick Zulian on 18/05/15.
//

#ifndef utopia_utopia_ASSIGN_HPP
#define utopia_utopia_ASSIGN_HPP

#include <string>
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Expression.hpp"
#include "utopia_Size.hpp"

namespace utopia {
    template<class Left, class Right>
    class Assign : public Expression< Assign<Left, Right> > {
    public:
        static const int Order = Left::Order;

        typedef typename Left::Scalar Scalar;
        
        Assign(Left &left, const Right &right) : _left(left), _right(right)
        {}

        inline Left &left() const { return _left; }
//        inline Left &left() const { return _left; }
        inline const Right &right() const { return _right; }

        std::string getClass() const {
            return "Assign<" + left().getClass() + ", " + right().getClass() + ">";
        }


    private:
        Left &_left;
        Right _right;
    };

    template<class Left, class Right>
    class Construct : public Expression< Construct<Left, Right> > {
    public:
        static const int Order = Left::Order;

        typedef typename Left::Scalar Scalar;

        Construct(Left &left, const Right &right) : _left(left), _right(right)
        {}

        inline Left &left() const { return _left; }
//        inline const Left &left() const { return _left; }
        inline const Right &right() const { return _right; }

        std::string getClass() const {
            return "Construct<" + left().getClass() + ", " + right().getClass() + ">";
        }

    private:
        Left &_left;
        Right _right;
    };

    template<class Left, class Right>
    inline Size size(const Assign<Left, Right> &expr)
    {
        return size(expr.right());
    }

    template<class Left, class Right>
    inline Size size(const Construct<Left, Right> &expr)
    {
        return size(expr.right());
    }

    template<class Left, class Right>
    inline Construct<Left, Right> construct(Expression<Left> &left, const Expression<Right> &right)
    {
        return Construct<Left, Right>(left.derived(), right.derived());
    }

    template<class Left, class Right>
    inline Assign<Left, Right> assign(Expression<Left> &left, const Expression<Right> &right)
    {
        return Assign<Left, Right>(left.derived(), right.derived());
    }



    template<class Left, class Right>
    class Traits< Construct<Left, Right> > : public Traits<Left> {};

    template<class Left, class Right>
    class Traits< Assign<Left, Right> >    : public Traits<Left> {};
}

#endif //utopia_utopia_ASSIGN_HPP
