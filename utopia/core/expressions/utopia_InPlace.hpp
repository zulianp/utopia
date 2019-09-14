//
// Created by Patrick Zulian on 22/05/15.
//

#ifndef UTOPIA_UTOPIA_INPLACE_HPP
#define UTOPIA_UTOPIA_INPLACE_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Expression.hpp"

namespace utopia {
    template<class Left, class Right, class Operation>
    class InPlace : public Expression< InPlace<Left, Right, Operation> > {
    public:
        InPlace(Left &left, const Right &right, const Operation operation = Operation()) : _left(left), _right(right), _operation(operation)
        {}

        inline Left &left() const { return _left; }
        inline const Right &right() const { return _right; }

        inline std::string get_class() const {
            return "InPlace<" + _left.get_class() + ", " + _right.get_class() + ", " + _operation.get_class() +  ">";
        }

        inline const Operation &operation() const { return _operation; }

    private:
        Left &_left;
        UTOPIA_STORE_CONST(Right) _right;
        Operation _operation;
    };


    template<class Left, class Right, class Operation>
    class Traits< InPlace<Left, Right, Operation> > : public Traits<Left> {};
}

#endif //UTOPIA_UTOPIA_INPLACE_HPP
