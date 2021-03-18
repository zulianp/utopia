#ifndef UTOPIA_LOCAL_TO_GLOBAL_HPP
#define UTOPIA_LOCAL_TO_GLOBAL_HPP

namespace utopia {
    template <class FunctionSpace, class ElementMatrices, class GlobalMatrix>
    class LocalToGlobal {};

    template <class FunctionSpace, class ElementMatrices, class GlobalMatrix>
    void local_to_global(const FunctionSpace &space, const ElementMatrices &element_matrices, GlobalMatrix &matrix) {
        LocalToGlobal<FunctionSpace, ElementMatrices, GlobalMatrix>::apply(space, element_matrices, matrix);
    }
}  // namespace utopia

#endif  // UTOPIA_LOCAL_TO_GLOBAL_HPP