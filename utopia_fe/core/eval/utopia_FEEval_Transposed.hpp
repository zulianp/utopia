#ifndef UTOPIA_FE_EVAL_TRANSPOSED_HPP
#define UTOPIA_FE_EVAL_TRANSPOSED_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_FEBackend.hpp"

namespace utopia {
    template<class Inner, class Traits, int Backend, int IsQuadData>
    class FEEval<Transposed<Inner>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::Transposed<Inner> Expr;

        inline static auto apply(
            const Expr &expr,
            AssemblyContext<Backend> &ctx) -> decltype(
                FEBackend<Backend>::transpose(FEEval<Inner, Traits, Backend, IsQuadData>::apply(expr.expr(), ctx), ctx)
            )
        {
            return FEBackend<Backend>::transpose(FEEval<Inner, Traits, Backend, IsQuadData>::apply(expr.expr(), ctx), ctx);
        }
    };

    template<class Expr, class AssemblyContext>
    class FunctionalTraits<Transposed<Expr>, AssemblyContext> {
    public:
        inline static int type(const Transposed<Expr> &expr, const AssemblyContext &ctx)
        {
            return FunctionalTraits<Expr, AssemblyContext>::type(expr.expr(), ctx);
        }

        inline static int order(const Transposed<Expr> &expr, const AssemblyContext &ctx)
        {
            return FunctionalTraits<Expr, AssemblyContext>::order(expr.expr(), ctx);
        }
    };
}

#endif //UTOPIA_FE_EVAL_TRANSPOSED_HPP
