//
// Created by Patrick Zulian on 30/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_UNARY_HPP
#define UTOPIA_UTOPIA_EVAL_UNARY_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {

    template<class Derived, typename Scalar, class Traits, int Backend>
    class Eval<Unary<
    Tensor<Derived, 1>,
    Reciprocal<Scalar> >,
    Traits, Backend> {
    public:
        inline static Derived apply(const Unary<Tensor<Derived, 1>, Reciprocal<Scalar> > &expr)
        {
            Derived result;

            UTOPIA_TRACE_BEGIN(expr);

            result.construct(Eval<Tensor<Derived, 1>, Traits>::apply(expr.expr()));
            result.transform(expr.operation());

            // //FIXME this is actually a binary thing
            // UTOPIA_BACKEND(Traits).apply_binary(
            //     result,
            //     expr.operation(),
            //     Eval<Tensor<Derived, 1>, Traits>::apply(expr.expr())
            // );

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template<class Expr, class Operation, class Traits, int Backend>
    class Eval<Unary<Expr, Operation>, Traits, Backend> {
    public:
        typedef typename TypeAndFill<Traits, Unary<Expr, Operation> >::Type Result;

        inline static Result apply(const Unary<Expr, Operation> &expr)
        {
            Result result;

            UTOPIA_TRACE_BEGIN(expr);

            result.construct(
                Eval<Expr, Traits>::apply(expr.expr())
            );

            result.transform(expr.operation());

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };
}

#endif //UTOPIA_UTOPIA_EVAL_UNARY_HPP
