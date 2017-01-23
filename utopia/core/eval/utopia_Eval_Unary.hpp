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
            //UTOPIA_EVENT_BEGIN(expr.getClass())
            Tensor result;
            const bool ok = UTOPIA_BACKEND(Traits).apply(
                    Eval<Wrapper<Tensor, 1>, Traits>::apply(expr.expr()),
                    expr.operation(),
                    result
            );

            assert(ok);

            //UTOPIA_EVENT_END(expr.getClass())
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
            const bool ok = UTOPIA_BACKEND(Traits).apply(
                    Eval<Expr, Traits>::apply(expr.expr()),
                    expr.operation(),
                    result
            );

            assert(ok);
            return result;
        }
    };
}

#endif //UTOPIA_UTOPIA_EVAL_UNARY_HPP
