#ifndef UTOPIA_FE_EVAL_WRAPPER_HPP
#define UTOPIA_FE_EVAL_WRAPPER_HPP

#include "utopia_FEForwardDeclarations.hpp"
#include "utopia_FunctionalTraits.hpp"

namespace utopia {
    template <class T, int Order, class AssemblyContext>
    class FunctionalTraits<Tensor<T, Order>, AssemblyContext> {
    public:
        inline static int type(const Tensor<T, Order> &expr, const AssemblyContext &) { return CONSTANT_FUNCTION; }

        inline static int order(const Tensor<T, Order> &expr, const AssemblyContext &) { return 0; }
    };
}  // namespace utopia

#endif  // UTOPIA_FE_EVAL_WRAPPER_HPP
