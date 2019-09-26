#ifndef UTOPIA_EVAL_NEGATE_HPP
#define UTOPIA_EVAL_NEGATE_HPP

#include "utopia_Eval_Unary.hpp"

namespace utopia {
    template<class Expr, class Traits, int Backend>
    class Eval<Negate<Expr>, Traits, Backend> {
    public:
        typedef utopia::Unary<Expr, Minus> DelegateT;
        typedef typename TypeAndFill<Traits, DelegateT>::Type Result;

        inline static Result apply(const Negate<Expr> &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);
            Result result = Eval<Expr, Traits>::apply(expr.expr());
            result.transform(Minus());
            UTOPIA_TRACE_END(expr);
            return result;
        }
    };
}

#endif //UTOPIA_EVAL_NEGATE_HPP
