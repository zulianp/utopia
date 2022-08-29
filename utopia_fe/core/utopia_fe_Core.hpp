#ifndef UTOPIA_FE_CORE_HPP
#define UTOPIA_FE_CORE_HPP

#include "utopia_FECoreForwardDeclarations.hpp"

#include "utopia_FEAssembler.hpp"

namespace utopia {
    template <class FunctionSpace>
    class Obstacle {};

    template <class FunctionSpace>
    class Contact {};

    template <class FunctionSpace>
    class FETransfer {};

    template <class FunctionSpace>
    class OmniAssembler {};

    template <class FunctionSpace>
    class IO {};

    enum AssemblyMode { ADD_MODE = 0, SUBTRACT_MODE = 1, OVERWRITE_MODE = 2 };

    template <class FE, class Op>
    class AssembleTraits {};

}  // namespace utopia

#endif  // UTOPIA_FE_CORE_HPP
