#ifndef UTOPIA_FE_EVAL_VALUE_FUNCTION_HPP
#define UTOPIA_FE_EVAL_VALUE_FUNCTION_HPP

#include "utopia_FEEval_Empty.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_FEBackend.hpp"
#include "utopia_FEFilter.hpp"
#include "utopia_FunctionalTraits.hpp"

namespace utopia {
    template<class In, class Out, int Order, class Fun, class Traits, int Backend, int IsQuadData>
    class FEEval<Filter<In, Out, Order, Fun>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::Filter<In, Out, Order, Fun> Expr;

        inline static auto apply(
            const Expr &expr,
            AssemblyContext<Backend> &ctx) -> decltype( 
            expr.eval(
                FEEval<In, Traits, Backend, IsQuadData>::apply(expr.expr(), ctx)
                ))
        {
            return expr.eval(
                FEEval<In, Traits, Backend, IsQuadData>::apply(expr.expr(), ctx)
            );
        }
    };

    template<class In, class Out, int Order, class Fun, class AssemblyContext>
    class FunctionalTraits<Filter<In, Out, Order, Fun>, AssemblyContext> {
    public:
        inline static int type(const Filter<In, Out, Order, Fun> &expr, const AssemblyContext &ctx)
        {
            return FunctionalTraits<In, AssemblyContext>::type(expr.expr(), ctx);
        }

        inline static int order(const Filter<In, Out, Order, Fun> &expr, const AssemblyContext &ctx)
        {
            return FunctionalTraits<In, AssemblyContext>::order(expr.expr(), ctx);
        }
    };
}


#endif //UTOPIA_FE_EVAL_VALUE_FUNCTION_HPP
