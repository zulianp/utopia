#ifndef UTOPIA_FE_EVAL_NORM2_HPP
#define UTOPIA_FE_EVAL_NORM2_HPP

#include "utopia_AssemblyContext.hpp"
#include "utopia_Eval_Empty.hpp"
#include "utopia_FEBackend.hpp"
#include "utopia_FEEval_Empty.hpp"
#include "utopia_Norm.hpp"

namespace utopia {

    template <class InnerExpr, int Type, class Traits, int Backend, int IsQuadData>
    class FEEval<Norm<InnerExpr, Type>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::Norm<InnerExpr, Type> Expr;

        inline static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx) -> decltype(
            FEBackend<Backend>::norm2(FEEval<InnerExpr, Traits, Backend, IsQuadData>::apply(expr.expr(), ctx), ctx)) {
            static_assert(Type == 2, "TODO other norms");

            return FEBackend<Backend>::norm2(FEEval<InnerExpr, Traits, Backend, IsQuadData>::apply(expr.expr(), ctx),
                                             ctx);
        }
    };

    template <class Inner, int Type, class AssemblyContext>
    class FunctionalTraits<Norm<Inner, Type>, AssemblyContext> {
    public:
        inline static int type(const Norm<Inner, Type> &expr, const AssemblyContext &ctx) {
            return FunctionalTraits<Inner, AssemblyContext>::type(expr.expr(), ctx);
        }

        inline static int order(const Norm<Inner, Type> &expr, const AssemblyContext &ctx) {
            return FunctionalTraits<Inner, AssemblyContext>::order(expr.expr(), ctx);
        }
    };
}  // namespace utopia

#endif  // UTOPIA_FE_EVAL_NORM2_HPP
