#ifndef UTOPIA_FE_EVAL_DETERMINANT_HPP
#define UTOPIA_FE_EVAL_DETERMINANT_HPP

#include "utopia_AssemblyContext.hpp"
#include "utopia_Eval_Empty.hpp"
#include "utopia_FEBackend.hpp"
#include "utopia_FindSpace.hpp"

namespace utopia {

    template <class InnerExpr, class Traits, int Backend, int IsQuadData>
    class FEEval<Determinant<InnerExpr>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::Determinant<InnerExpr> Expr;

        inline static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx) -> decltype(
            FEBackend<Backend>::determinant(FEEval<InnerExpr, Traits, Backend, IsQuadData>::apply(expr.expr(), ctx),
                                            ctx)) {
            return FEBackend<Backend>::determinant(
                FEEval<InnerExpr, Traits, Backend, IsQuadData>::apply(expr.expr(), ctx), ctx);
        }
    };

    template <class Inner, class AssemblyContext>
    class FunctionalTraits<Determinant<Inner>, AssemblyContext> {
    public:
        inline static int type(const Determinant<Inner> &expr, const AssemblyContext &ctx) {
            return FunctionalTraits<Inner, AssemblyContext>::type(expr.expr(), ctx);
        }

        inline static int order(const Determinant<Inner> &expr, const AssemblyContext &ctx) {
            return FunctionalTraits<Inner, AssemblyContext>::order(expr.expr(), ctx);
        }
    };
}  // namespace utopia

#endif  // UTOPIA_FE_EVAL_DETERMINANT_HPP
