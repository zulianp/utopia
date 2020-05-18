//
// Created by Patrick Zulian on 28/05/15.
//

#ifndef UTOPIA_UTOPIA_MULTIPLY_HPP
#define UTOPIA_UTOPIA_MULTIPLY_HPP

#include "utopia_Castable.hpp"
#include "utopia_Expression.hpp"
#include "utopia_Size.hpp"

// #define TENSOR_ORDER_MULTIPLY(Left_, Right_) (Left_::Order < Right_::Order)? Left_::Order : Right_::Order

namespace utopia {
    template <class Left, class Right>
    class MultiplyTensorOrder {
    public:
        static const int Order =
            (Left::Order == 0)
                ? Right::Order
                : ((Right::Order == 0) ? Left::Order : (Left::Order < Right::Order) ? Left::Order : Right::Order);
    };

#define TENSOR_ORDER_MULTIPLY(Left_, Right_) (MultiplyTensorOrder<Left_, Right_>::Order)

    template <class _Left, class _Right>
    class Multiply : public Expression<Multiply<_Left, _Right> >,
                     public Castable<Multiply<_Left, _Right>, TENSOR_ORDER_MULTIPLY(_Left, _Right)> {
    public:
        using Left = _Left;
        using Right = _Right;
        using Scalar = decltype(typename Traits<Left>::Scalar() * typename Traits<Right>::Scalar());

        static const int Order = TENSOR_ORDER_MULTIPLY(_Left, _Right);

        Multiply(const Left &left, const Right &right) : _left(left), _right(right) {}

        const Left &left() const { return _left; }
        const Right &right() const { return _right; }

        std::string get_class() const override {
            return "Multiply<" + _left.get_class() + ", " + _right.get_class() + ">";
        }

        virtual ~Multiply() = default;

    private:
        UTOPIA_STORE_CONST(Left) _left;
        UTOPIA_STORE_CONST(Right) _right;
    };

    template <class Left, class Right>
    inline Size size(const Multiply<Left, Right> &expr) {
        Size result(2);
        result.set(0, size(expr.left()).get(0));

        Size r_size = size(expr.right());
        if (r_size.n_dims() == 1) {
            result.set(1, 1);
        } else {
            result.set(1, size(expr.right()).get(1));
        }

        return result;
    }

    template <class Left, class Right>
    class Traits<Multiply<Left, Right> > : public Traits<typename ChooseType<Left, Right, Right>::Type> {};
}  // namespace utopia

#undef TENSOR_ORDER_MULTIPLY
#endif  // UTOPIA_UTOPIA_MULTIPLY_HPP
