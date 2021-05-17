#ifndef UTOPIA_FE_EVAL_FACTORY_HPP
#define UTOPIA_FE_EVAL_FACTORY_HPP

#include "utopia_AssemblyContext.hpp"
#include "utopia_Eval_Empty.hpp"
#include "utopia_Factory.hpp"

namespace utopia {

    template <class Type, int Order, class Traits, int Backend, int IsQuadData>
    class FEEval<Factory<Type, Order>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::Factory<Type, Order> Expr;

        inline static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx)
            -> decltype(FEBackend<Backend>::build(expr.size(), expr.type(), ctx)) {
            return FEBackend<Backend>::build(expr.size(), expr.type(), ctx);
        }
    };

    template <class Type, int Order, class AssemblyContext>
    class FunctionalTraits<Factory<Type, Order>, AssemblyContext> {
    public:
        inline static int type(const Factory<Type, Order> &expr, const AssemblyContext &ctx) {
            return CONSTANT_FUNCTION;
        }

        inline static int order(const Factory<Type, Order> &expr, const AssemblyContext &ctx) { return 0; }
    };
}  // namespace utopia

#endif  // UTOPIA_FE_EVAL_FACTORY_HPP
