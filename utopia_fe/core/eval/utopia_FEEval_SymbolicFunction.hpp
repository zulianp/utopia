#ifndef UTOPIA_FE_EVAL_SYMBOLIC_FUNCTION_HPP
#define UTOPIA_FE_EVAL_SYMBOLIC_FUNCTION_HPP

#ifdef UTOPIA_WITH_TINY_EXPR

#include "utopia_AssemblyContext.hpp"
#include "utopia_Eval_Empty.hpp"
#include "utopia_FEBackend.hpp"
#include "utopia_SymbolicFunction.hpp"

namespace utopia {
    template <class Traits, int Backend, int IsQuadData>
    class FEEval<SymbolicFunction, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::SymbolicFunction Expr;

        inline static auto apply(const SymbolicFunction &expr, AssemblyContext<Backend> &ctx)
            -> decltype(FEBackend<Backend>::apply(expr, ctx)) {
            return FEBackend<Backend>::apply(expr, ctx);
        }
    };

    template <class AssemblyContext>
    class FunctionalTraits<SymbolicFunction, AssemblyContext> {
    public:
        inline static int type(const SymbolicFunction &expr, const AssemblyContext &ctx) { return 0; }

        inline static int order(const SymbolicFunction &expr, const AssemblyContext &ctx) { return 0; }
    };
}  // namespace utopia

#endif  // UTOPIA_WITH_TINY_EXPR
#endif  // UTOPIA_FE_EVAL_SYMBOLIC_FUNCTION_HPP
