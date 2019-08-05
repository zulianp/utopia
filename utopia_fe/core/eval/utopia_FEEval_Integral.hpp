#ifndef UTOPIA_FE_EVAL_INTEGRAL_HPP
#define UTOPIA_FE_EVAL_INTEGRAL_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_FEBackend.hpp"

namespace utopia {
    template<class Inner, class Traits, int Backend, int IsQuadData>
    class FEEval< Integral<Inner>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::Integral<Inner> Expr;

        inline static auto apply(
            const Expr &expr,
            AssemblyContext<Backend> &ctx) -> decltype(
                FEBackend<Traits::Backend>::integrate(
                                FEEval<Inner, Traits, Backend, IsQuadData>::apply(expr.expr(), ctx),
                            ctx)
            )
        {
            //Check subtree properties
            //TODO
            return FEBackend<Traits::Backend>::integrate(
                FEEval<Inner, Traits, Backend, IsQuadData>::apply(expr.expr(), ctx),
            ctx);
        }
    };
}

#endif //UTOPIA_FE_EVAL_INTEGRAL_HPP
