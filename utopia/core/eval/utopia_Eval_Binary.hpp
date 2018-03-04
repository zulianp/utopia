//
// Created by Patrick Zulian on 29/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_BINARY_HPP_HPP
#define UTOPIA_UTOPIA_EVAL_BINARY_HPP_HPP

#include "utopia_Eval_Empty.hpp"


namespace utopia {

    template<class ScalarT, class Right, class Operation, class Traits, int Backend>
    class Eval<Binary<Number<ScalarT>, Right, Operation>, Traits, Backend> {
    public:
        typedef typename TypeAndFill<Traits, Binary<Number<ScalarT>, Right, Operation> >::Type Result;
        inline static Result apply(const Binary<Number<ScalarT>, Right, Operation> &expr)
        {
            Result result;

            UTOPIA_LOG_BEGIN(expr);

            UTOPIA_BACKEND(Traits).apply_binary(
                result,
                expr.left(),
                expr.operation(),
                Eval<Right, Traits>::apply(expr.right())
                );

            UTOPIA_LOG_END(expr);
            return result;
        }
    };

    template<class Left, class Right, class Operation, class Traits, int Backend>
    class Eval<Binary<Number<Left>, Number<Right>, Operation>, Traits, Backend> {
    public:

        inline static auto apply(const Binary<Number<Left>, Number<Right>, Operation> &expr) -> decltype(Left() + Right())
        {
            Left l = expr.left();
            Right r = expr.right();
            return expr.operation().apply(l, r);
        }
    };

    template<class Left, class Right, class Operation, class Traits, int Backend>
    class Eval<Binary<Left, Right, Operation>, Traits, Backend> {
    public:
        typedef typename utopia::TypeAndFill<Traits, Binary<Left, Right, Operation> >::Type Result;

        inline static Result apply(const Binary<Left, Right, Operation> &expr) {
            Result result;

            UTOPIA_LOG_BEGIN(expr);

            UTOPIA_BACKEND(Traits).apply_binary(
                result,
                Eval<Left,  Traits>::apply(expr.left()),
                expr.operation(),
                Eval<Right, Traits>::apply(expr.right())                 
                );

            UTOPIA_LOG_END(expr);
            return result;
        }
    };

    template<class Left, class ScalarT, class Traits, int Backend>
    class Eval<Binary<Left, Number<ScalarT>, Multiplies>, Traits, Backend> {
    public:

        typedef typename TypeAndFill<Traits, Binary<Left, Number<ScalarT>, Multiplies> >::Type Result;
        inline static Result apply(const Binary<Left, Number<ScalarT>, Multiplies> &expr)
        {
            Result result;

            UTOPIA_LOG_BEGIN(expr);

            UTOPIA_BACKEND(Traits).apply_binary(
                result,
                expr.right(),
                expr.operation(),
                Eval<Left, Traits>::apply(expr.left())
                );

            UTOPIA_LOG_END(expr);
            return result;
        }
    };


    template<class Left, class Right, class Traits, int Backend>
    class Eval<OuterProduct<Left, Right>, Traits, Backend> {
    public:
        typedef utopia::OuterProduct<Left, Right> Expr;

        inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr) {
            EXPR_TYPE(Traits, Expr) result;

            UTOPIA_LOG_BEGIN(expr);

            UTOPIA_BACKEND(Traits).kronecker_product(
                result,
                Eval<Left, Traits>::apply(expr.left()),
                Eval<Right, Traits>::apply(expr.right())
                );

            UTOPIA_LOG_END(expr);
            return result;
        }
    };

}

#endif //UTOPIA_UTOPIA_EVAL_BINARY_HPP_HPP
