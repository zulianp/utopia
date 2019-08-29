#ifndef UTOPIA_FE_EVAL_DIV_HPP
#define UTOPIA_FE_EVAL_DIV_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_FEBackend.hpp"

namespace utopia {

    template<class Tensor, class Traits, int Backend, int IsQuadData>
    class FEEval<Divergence<Tensor>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::Divergence<Tensor> Expr;

        inline static auto apply(
            const Divergence<Tensor> &expr,
            AssemblyContext<Backend> &ctx) -> decltype( FEBackend<Backend>::div(expr.expr(), ctx) )
        {
            return FEBackend<Backend>::div(expr.expr(), ctx);
        }
    };

}

#endif //UTOPIA_FE_EVAL_DIV_HPP
