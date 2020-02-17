#ifndef UTOPIA_EVAL_NEGATE_HPP
#define UTOPIA_EVAL_NEGATE_HPP

#include "utopia_Eval_Unary.hpp"

namespace utopia {
    template<class InnerExpr, class Traits, int Backend>
    class Eval<Negate<InnerExpr>, Traits, Backend> {
    public:
        using Expr = utopia::Negate<InnerExpr>;
        using Result = EXPR_TYPE(Traits, Expr);

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Negate<InnerExpr> &expr, Result &result)
        {
            UTOPIA_TRACE_BEGIN(expr);
            result = Eval<InnerExpr, Traits>::apply(expr.expr());
            result.transform(Minus());
            UTOPIA_TRACE_END(expr);
        }
    };
}

#endif //UTOPIA_EVAL_NEGATE_HPP
