#ifndef utopia_utopia_BINARY_HPP
#define utopia_utopia_BINARY_HPP

#include "utopia_Castable.hpp"
#include "utopia_Expression.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Literal.hpp"
#include "utopia_Size.hpp"
#include "utopia_Traits.hpp"

#define TENSOR_ORDER_BINARY(Left_, Right_) (Left_::Order > Right_::Order) ? Left_::Order : Right_::Order

namespace utopia {

    template <class _Left, class _Right, class _Operation>
    class Binary : public Expression<Binary<_Left, _Right, _Operation> >,
                   public Castable<Binary<_Left, _Right, _Operation>, TENSOR_ORDER_BINARY(_Left, _Right)> {
    public:
        using Left = _Left;
        using Right = _Right;
        using Operation = _Operation;
        using Scalar = decltype(typename Left::Scalar() + typename Right::Scalar());

        static const int Order = TENSOR_ORDER_BINARY(_Left, _Right);

        Binary(const Left &left, const Right &right, const Operation operation = Operation())
            : _left(left), _right(right), _operation(operation) {}

        const Left &left() const { return _left; }
        const Right &right() const { return _right; }
        const Operation &operation() const { return _operation; }

        std::string get_class() const override {
            return GetClass<Operation>() + "<" + _left.get_class() + ", " + _right.get_class() + ">";
        }

        ~Binary() = default;

    private:
        UTOPIA_STORE_CONST(Left) _left;
        UTOPIA_STORE_CONST(Right) _right;
        Operation _operation;
    };

    template <class Left, class Right, class Operation>
    class Traits<Binary<Left, Right, Operation> > : public Traits<typename MostDescriptive<Left, Right>::Type> {};

    template <class Left, class Right, class Operation>
    Size size(const Binary<Left, Right, Operation> &expr) {
        if (Left::Order > Right::Order) {
            return size(expr.left());
        } else {
            return size(expr.right());
        }
    }

}  // namespace utopia

#undef TENSOR_ORDER_BINARY

#endif  // utopia_utopia_BINARY_HPP
