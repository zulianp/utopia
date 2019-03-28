#ifndef UTOPIA_FE_EVAL_CONTEXT_FUNCTION_HPP
#define UTOPIA_FE_EVAL_CONTEXT_FUNCTION_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_FEBackend.hpp"


namespace utopia {
    template<class Out, class Fun, class Traits, int Backend, int IsQuadData>
    class FEEval<ContextFunction<Out, Fun>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::ContextFunction<Out, Fun> Expr;

        inline static auto apply(
            const Expr &expr,
            AssemblyContext<Backend> &ctx) -> Out
        {
            return expr.eval(ctx);
        }
    };



    template<class Out, class Fun, class AssemblyContext>
    class FunctionalTraits<ContextFunction<Out, Fun>, AssemblyContext> {
    public:
        inline static int type(const ContextFunction<Out, Fun> &expr, const AssemblyContext &ctx)
        {
            return 0;
        }

        inline static int order(const ContextFunction<Out, Fun> &expr, const AssemblyContext &ctx)
        {
            return 0;
        }
    };
}


#endif //UTOPIA_FE_EVAL_CONTEXT_FUNCTION_HPP
