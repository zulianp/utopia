#ifndef UTOPIA_EVAL_TENSOR_REDUCE_HPP
#define UTOPIA_EVAL_TENSOR_REDUCE_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Operators.hpp"

namespace utopia {

    template<class Expr, class Operation, class Traits, int Backend>
    class Eval< TensorReduce<Expr, Operation>, Traits, Backend> {
    public:
        typedef typename TypeAndFill<Traits,TensorReduce<Expr, Operation> >::Type Result;

        inline static Result apply(const TensorReduce<Expr, Operation> &expr) {
            Result result;

            UTOPIA_TRACE_BEGIN(expr);

            // UTOPIA_BACKEND(Traits).apply_tensor_reduce(
            //         result,
            //         Eval<Expr,  Traits>::apply(expr.expr()),
            //         expr.operation(),
            //         expr.dim()
            // );

            apply_aux(
                Eval<Expr,  Traits>::apply(expr.expr()),
                expr.operation(),
                expr.dim(),
                result
            );

            UTOPIA_TRACE_END(expr);
            return result;
        }

        template<class Operand>
        inline static void apply_aux(const Operand &t, const Plus &, const int dim, Result &result)
        {
            if(dim == 1) {
                t.row_sum(result);
            } else {
                t.col_sum(result);
            }
        }

        template<class Operand>
        inline static void apply_aux(const Operand &t, const Max &, const int dim, Result &result)
        {
            if(dim == 1) {
                t.row_max(result);
            } else {
                t.col_max(result);
            }
        }

        template<class Operand>
        inline static void apply_aux(const Operand &t, const Min &, const int dim, Result &result)
        {
            if(dim == 1) {
                t.row_min(result);
            } else {
                t.col_min(result);
            }
        }
    };
}


#endif //UTOPIA_EVAL_TENSOR_REDUCE_HPP
