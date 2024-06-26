#ifndef UTOPIA_EVAL_RESIDUAL_HPP
#define UTOPIA_EVAL_RESIDUAL_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_ForwardDeclarations.hpp"

namespace utopia {
    template <class A, class X, class B>
    using ResidualExpr = Binary<Tensor<B, 1>, Multiply<Tensor<A, 2>, Tensor<X, 1>>, Minus>;

    template <class A, class X, class B, class Traits, int Backend>
    class Eval<ResidualExpr<A, X, B>, Traits, Backend> {
    public:
        typedef utopia::ResidualExpr<A, X, B> Expr;
        using Result = X;

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result) {
            UTOPIA_TRACE_BEGIN(expr);

            const auto &a = expr.right().left().derived();
            const auto &x = expr.right().right().derived();
            const auto &b = expr.left().derived();

            if (result.is_alias(x) || result.is_alias(b)) {
                X temp;
                a.multiply(x, temp);

                if (result.is_alias(b)) {
                    result.axpy(-1.0, temp);
                } else {
                    result = b;
                    result.axpy(-1.0, temp);
                }
            } else {
                a.multiply(x, result);
                result.scale(-1.0);
                result.axpy(1.0, b);
            }

            UTOPIA_TRACE_END();
        }
    };

    // template<class Left, class A, class X, class B, class Traits, int Backend>
    // class Eval<Assign<Left, ResidualExpr<A, X, B>>, Traits, Backend> {
    // public:
    //     typedef utopia::ResidualExpr<A, X, B> Right;
    //     typedef utopia::Assign<Left, Right> Expr;

    //     inline static void apply(const Expr &assign_expr) {
    //         UTOPIA_TRACE_BEGIN(assign_expr);
    //         auto &&expr = assign_expr.right();

    //         auto &r = assign_expr.left().derived();
    //         const auto &a = expr.right().left().derived();
    //         const auto &x = expr.right().right().derived();
    //         const auto &b = expr.left().derived();

    //         if(r.is_alias(x) || r.is_alias(b)) {
    //             X temp;
    //             a.multiply(x, temp);

    //             if(r.is_alias(b)) {
    //                 r.axpy(-1.0, temp);
    //             } else {
    //                 r = b;
    //                 r.axpy(-1.0, temp);
    //             }
    //         } else {
    //             a.multiply(x, r);
    //             r.scale(-1.0);
    //             r.axpy(1.0, b);
    //         }

    //         UTOPIA_TRACE_END(assign_expr);
    //     }
    // };
}  // namespace utopia

#endif  // UTOPIA_EVAL_RESIDUAL_HPP
