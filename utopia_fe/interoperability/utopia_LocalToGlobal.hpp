#ifndef UTOPIA_LOCAL_TO_GLOBAL_HPP
#define UTOPIA_LOCAL_TO_GLOBAL_HPP

#include "utopia_fe_Core.hpp"

namespace utopia {
    template <class FunctionSpace, class ElementTensors, class GlobalTensor, typename... Args>
    class LocalToGlobal {};

    template <class FunctionSpace, class ElementTensors, class GlobalTensor, typename... Args>
    void local_to_global(const FunctionSpace &space,
                         const ElementTensors &element_matrices,
                         AssemblyMode mode,
                         GlobalTensor &tensor,
                         Args &&... args) {
        LocalToGlobal<FunctionSpace, ElementTensors, GlobalTensor>::apply(
            space, element_matrices, mode, tensor, args...);
    }

    template <class FunctionSpace, class GlobalTensor, class ElementTensors, typename... Args>
    class GlobalToLocal {};

    template <class FunctionSpace, class GlobalTensor, class ElementTensors, typename... Args>
    void global_to_local(const FunctionSpace &space,
                         const GlobalTensor &tensor,
                         AssemblyMode mode,
                         ElementTensors &element_matrices,
                         Args &&... args) {
        GlobalToLocal<FunctionSpace, GlobalTensor, ElementTensors>::apply(
            space, tensor, mode, element_matrices, args...);
    }
}  // namespace utopia

#endif  // UTOPIA_LOCAL_TO_GLOBAL_HPP