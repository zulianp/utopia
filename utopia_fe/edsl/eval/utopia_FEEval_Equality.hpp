#ifndef UTOPIA_FE_EVAL_EQUALITY_HPP
#define UTOPIA_FE_EVAL_EQUALITY_HPP

#include "utopia_FEForwardDeclarations.hpp"

namespace utopia {
    template <class Left, class Right, class AssemblyContext>
    class FunctionalTraits<Equality<Left, Right>, AssemblyContext> {
    public:
        inline static int type(const Equality<Left, Right> &expr, const AssemblyContext &ctx) {
            return std::max(FunctionalTraits<Left, AssemblyContext>::type(expr.left(), ctx),
                            FunctionalTraits<Right, AssemblyContext>::type(expr.right(), ctx));
        }

        inline static int order(const Equality<Left, Right> &expr, const AssemblyContext &ctx) {
            return std::max(FunctionalTraits<Left, AssemblyContext>::order(expr.left(), ctx),
                            FunctionalTraits<Right, AssemblyContext>::order(expr.right(), ctx));
        }
    };
}  // namespace utopia

#endif  // UTOPIA_FE_EVAL_EQUALITY_HPP
