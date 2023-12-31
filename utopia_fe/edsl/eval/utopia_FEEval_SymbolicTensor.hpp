#ifndef UTOPIA_FE_EVAL_SYMBOLIC_TENSOR_HPP
#define UTOPIA_FE_EVAL_SYMBOLIC_TENSOR_HPP

#include "utopia_AssemblyContext.hpp"
#include "utopia_Eval_Empty.hpp"
#include "utopia_Factory.hpp"

namespace utopia {

    // FIXME
    template <class Left, class Type, int Order, class Traits, int Backend, int IsQuadData>
    class FEEval<Multiply<Trace<Left>, SymbolicTensor<Type, Order> >, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::Multiply<Trace<Left>, SymbolicTensor<Type, Order> > Expr;

        inline static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx)
            -> decltype(FEBackend<Backend>::trace_times_identity(
                FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left().expr(), ctx),
                ctx)

            ) {
            return FEBackend<Backend>::trace_times_identity(
                FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left().expr(), ctx), ctx);
        }
    };

    // FIXME
    template <class Type, int Order, class Traits, int Backend, int IsQuadData>
    class FEEval<SymbolicTensor<Type, Order>, Traits, Backend, IsQuadData> {
    public:
        inline static SymbolicTensor<Type, Order> &apply(const SymbolicTensor<Type, Order> &expr,
                                                         AssemblyContext<Backend> &ctx) {
            return expr;
        }
    };

    template <class Type, int Order, class AssemblyContext>
    class FunctionalTraits<SymbolicTensor<Type, Order>, AssemblyContext> {
    public:
        inline static int type(const SymbolicTensor<Type, Order> &expr, const AssemblyContext &ctx) {
            return CONSTANT_FUNCTION;
        }

        inline static int order(const SymbolicTensor<Type, Order> &expr, const AssemblyContext &ctx) { return 0; }
    };
}  // namespace utopia

#endif  // UTOPIA_FE_EVAL_SYMBOLIC_TENSOR_HPP
