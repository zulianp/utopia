#ifndef UTOPIA_UTOPIA_EVAL_TRANSPOSE_HPP
#define UTOPIA_UTOPIA_EVAL_TRANSPOSE_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {

    template<class Tensor, class Traits, int Backend>
    class Eval<Transposed<Tensor>, Traits, Backend> {
    public:
        using Expr = utopia::Transposed<Tensor>;
        using Result = EXPR_TYPE(Traits, Expr);
        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &t, Result &result)
        {
            UTOPIA_TRACE_BEGIN(t);
            Eval<Tensor, Traits>::apply(t.expr()).transpose(result);
            UTOPIA_TRACE_END(t);
        }
    };
}

#endif //UTOPIA_UTOPIA_EVAL_TRANSPOSE_HPP
