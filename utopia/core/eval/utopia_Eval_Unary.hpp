//
// Created by Patrick Zulian on 30/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_UNARY_HPP
#define UTOPIA_UTOPIA_EVAL_UNARY_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {

    template <class Derived, typename Scalar, class Traits, int Backend>
    class Eval<Unary<Tensor<Derived, 1>, Reciprocal<Scalar> >, Traits, Backend> {
    public:
        using Expr = utopia::Unary<Tensor<Derived, 1>, Reciprocal<Scalar> >;
        using Result = Derived;

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result) {
            UTOPIA_TRACE_BEGIN(expr);
            auto &&t = Eval<Tensor<Derived, 1>, Traits>::apply(expr.expr());

            if (result.is_alias(t)) {
                result.transform(expr.operation());
            } else {
                result.construct(t);
                result.transform(expr.operation());
            }

            UTOPIA_TRACE_END(expr);
        }
    };

    template <class InnerExpr, class Operation, class Traits, int Backend>
    class Eval<Unary<InnerExpr, Operation>, Traits, Backend> {
    public:
        using Expr = utopia::Unary<InnerExpr, Operation>;
        using Result = EXPR_TYPE(Traits, Expr);

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result) {
            UTOPIA_TRACE_BEGIN(expr);
            auto &&t = Eval<InnerExpr, Traits>::apply(expr.expr());

            if (result.is_alias(t)) {
                result.transform(expr.operation());
            } else {
                result.construct(t);
                result.transform(expr.operation());
            }

            UTOPIA_TRACE_END(expr);
        }
    };
}  // namespace utopia

#endif  // UTOPIA_UTOPIA_EVAL_UNARY_HPP
