#ifndef UTOPIA_EVAL_TENSOR_REDUCE_HPP
#define UTOPIA_EVAL_TENSOR_REDUCE_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Operators.hpp"

namespace utopia {

    template <class InnerExpr, class Operation, class Traits, int Backend>
    class Eval<TensorReduce<InnerExpr, Operation>, Traits, Backend> {
    public:
        using Expr = utopia::TensorReduce<InnerExpr, Operation>;
        using Result = EXPR_TYPE(Traits, Expr);

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result) {
            UTOPIA_TRACE_BEGIN(expr);

            apply_aux(Eval<InnerExpr, Traits>::apply(expr.expr()), expr.operation(), expr.dim(), result);

            UTOPIA_TRACE_END(expr);
        }

        template <class Operand>
        inline static void apply_aux(const Operand &t, const Plus &, const int dim, Result &result) {
            if (dim == 1) {
                t.row_sum(result);
            } else {
                t.col_sum(result);
            }
        }

        template <class Operand>
        inline static void apply_aux(const Operand &t, const Max &, const int dim, Result &result) {
            if (dim == 1) {
                t.row_max(result);
            } else {
                t.col_max(result);
            }
        }

        template <class Operand>
        inline static void apply_aux(const Operand &t, const Min &, const int dim, Result &result) {
            if (dim == 1) {
                t.row_min(result);
            } else {
                t.col_min(result);
            }
        }
    };
}  // namespace utopia

#endif  // UTOPIA_EVAL_TENSOR_REDUCE_HPP
