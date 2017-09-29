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

            const bool ok =UTOPIA_BACKEND(Traits).apply(
                    expr.left(),
                    Eval<Right, Traits>::apply(expr.right()),
                    expr.operation(),
                    result);

            ASSERT(ok);


            // [cleaup] This eval maybe can be removed completely
            // [new backend map concept]
            // [optimized][minimal] backend
            // UTOPIA_BACKEND(Traits).apply(result, alpha, Operation, right);

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

            UTOPIA_BACKEND(Traits).outer(
                    Eval<Left, Traits>::apply(expr.left()),
                    Eval<Right, Traits>::apply(expr.right()),
                    result
            );

            // [new backend map concept]
            // [optimized][minimal] backend
            // UTOPIA_BACKEND(Traits).apply(result, left, Kron, right);

            UTOPIA_LOG_END(expr);
            return result;
        }
    };


    template<class Left, class Right, class Operation, class Traits, int Backend>
    class Eval<Binary<Left, Right, Operation>, Traits, Backend> {
    public:
        typedef typename utopia::TypeAndFill<Traits, Binary<Left, Right, Operation> >::Type Result;

        inline static Result apply(const Binary<Left, Right, Operation> &expr) {
            Result result;

            UTOPIA_LOG_BEGIN(expr);

            const bool ok = UTOPIA_BACKEND(Traits).apply(
                    Eval<Left,  Traits>::apply(expr.left()),
                    Eval<Right, Traits>::apply(expr.right()),
                    expr.operation(),
                    result
            );

            ASSERT(ok);

            // [new backend map concept]
            // [optimized][minimal] backend
            // UTOPIA_BACKEND(Traits).apply(result, left, Operation, right);

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

            bool ok = UTOPIA_BACKEND(Traits).apply(
                    expr.right(),
                    Eval<Left, Traits>::apply(expr.left()),
                    expr.operation(),
                    result);

            ASSERT(ok);

            // [cleaup] This eval maybe can be removed completely
            // [new backend map concept]
            // [optimized][minimal] backend
            // UTOPIA_BACKEND(Traits).apply(result, alpha, Multiples, right);

			UTOPIA_LOG_END(expr);
            return result;
        }
    };
}

#endif //UTOPIA_UTOPIA_EVAL_BINARY_HPP_HPP
