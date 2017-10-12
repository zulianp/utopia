//
// Created by Patrick Zulian on 15/05/15.
//

#ifndef utopia_utopia_BINARY_HPP
#define utopia_utopia_BINARY_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Expression.hpp"
#include "utopia_Size.hpp"
#include "utopia_Literal.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Castable.hpp"

#define TENSOR_ORDER_BINARY(Left_, Right_) (Left_::Order > Right_::Order)? Left_::Order : Right_::Order

namespace utopia {

    template<class _Left, class _Right, class _Operation>
    class Binary : public Expression<Binary<_Left, _Right, _Operation> >, 
                   public Castable<Binary<_Left, _Right, _Operation>,
                                   TENSOR_ORDER_BINARY(_Left, _Right)> {
    public:
        typedef _Left Left;
        typedef _Right Right;
        typedef _Operation Operation;
        typedef decltype(typename Left::Scalar() + typename Right::Scalar()) Scalar;

        static const int Order = TENSOR_ORDER_BINARY(_Left, _Right);

        Binary(const Left &left, const Right &right, const Operation operation = Operation())
                : _left(left), _right(right), _operation(operation)
        {}

        const Left &left() const { return _left; }
        const Right &right() const { return _right; }
        const Operation &operation() const { return _operation; }

        virtual std::string getClass() const
        {
            return GetClass<Operation>() + "<" + _left.getClass() + ", " + _right.getClass() + ">";
        }

        virtual ~Binary() { }

    private:
        UTOPIA_STORE_CONST(Left) _left;
        UTOPIA_STORE_CONST(Right) _right;
        Operation _operation;
    };

    template<class Left, class Right, class Operation>
    class Traits< Binary<Left, Right, Operation> > : public Traits< typename MostDescriptive<Left, Right>::Type > {};

    template<class Left, class Right, class Operation>
    Size size(const Binary<Left, Right, Operation> &expr)
    {
        if(Left::Order > Right::Order) {
            return size(expr.left());
        } else {
            return size(expr.right());
        }
    }
}

#undef TENSOR_ORDER_BINARY

#endif //utopia_utopia_BINARY_HPP
