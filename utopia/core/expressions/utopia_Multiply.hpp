//
// Created by Patrick Zulian on 28/05/15.
//

#ifndef UTOPIA_UTOPIA_MULTIPLY_HPP
#define UTOPIA_UTOPIA_MULTIPLY_HPP

#include "utopia_Expression.hpp"
#include "utopia_Size.hpp"
#include "utopia_Castable.hpp"

#define TENSOR_ORDER_MULTIPLY(Left_, Right_) (static_cast<long>(Left_::Order) < static_cast<long>(Right_::Order))? \
		static_cast<long>(Left_::Order) : static_cast<long>(Right_::Order)

namespace utopia {
    template<class _Left, class _Right>
    class Multiply : public Expression< Multiply<_Left, _Right> >, 
                     public Castable< Multiply<_Left, _Right>,
                                      TENSOR_ORDER_MULTIPLY(_Left, _Right)>
                    {
    public:
        typedef _Left Left;
        typedef _Right Right;
        typedef decltype(typename Left::Scalar() * typename Right::Scalar()) Scalar;

        enum {
            Order = TENSOR_ORDER_MULTIPLY(_Left, _Right)
        };

        Multiply(const Left &left, const Right &right)
                : _left(left), _right(right)
        {}

        const Left &left() const { return _left; }
        const Right &right() const { return _right; }

        virtual std::string getClass() const
        {
            return  "Multiply<" + _left.getClass() + ", " + _right.getClass() + ">";
        }

        virtual ~Multiply() { }

    private:
        UTOPIA_STORE_CONST(Left) _left;
        UTOPIA_STORE_CONST(Right) _right;
    };

    template<class Left, class Right>
    inline Size size(const Multiply<Left, Right> &expr)
    {
        Size result(2);
        result.set(0, size(expr.left()).get(0));

        Size r_size = size(expr.right());
        if(r_size.n_dims() == 1) {
            result.set(1, 1);
        } else {
            result.set(1, size(expr.right()).get(1));
        }
        
        return result;
    }

    template<class Left, class Right>
    class Traits< Multiply<Left, Right> > : public Traits<typename ChooseType<Left, Right, Right>::Type > {};
}

#undef TENSOR_ORDER_MULTIPLY
#endif //UTOPIA_UTOPIA_MULTIPLY_HPP
