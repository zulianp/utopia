//
// Created by Patrick Zulian on 30/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_UNARY_HPP
#define UTOPIA_UTOPIA_EVAL_UNARY_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {

    template<class Tensor, typename Scalar, class Traits, int Backend>
    class Eval<Unary<
            Wrapper<Tensor, 1>,
            Reciprocal<Scalar> >,
            Traits, Backend> {
    public:

        inline static Tensor apply(const Unary<Wrapper<Tensor, 1>, Reciprocal<Scalar> > &expr)
        {
            Tensor result;

            UTOPIA_LOG_BEGIN(expr);

            const bool ok = UTOPIA_BACKEND(Traits).apply(
                    Eval<Wrapper<Tensor, 1>, Traits>::apply(expr.expr()),
                    expr.operation(),
                    result
            );

            ASSERT(ok);

            UTOPIA_LOG_END(expr);
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

            UTOPIA_LOG_BEGIN(expr);

            const bool ok = UTOPIA_BACKEND(Traits).apply(
                    Eval<Expr, Traits>::apply(expr.expr()),
                    expr.operation(),
                    result
            );

            ASSERT(ok);

            UTOPIA_LOG_END(expr);
            return result;
        }
    };
}

#endif //UTOPIA_UTOPIA_EVAL_UNARY_HPP
