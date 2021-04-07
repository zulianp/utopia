#ifndef UTOPIA_FE_CORE_HPP
#define UTOPIA_FE_CORE_HPP

namespace utopia {
    template <class FunctionSpace>
    class Obstacle {};

    template <class FunctionSpace>
    class FETransfer {};

    template <class FunctionSpace>
    class OmniAssembler {};

    enum AssemblyMode { ADD_MODE = 0, SUBTRACT_MODE = 1, OVERWRITE_MODE = 2 };
}  // namespace utopia

#endif  // UTOPIA_FE_CORE_HPP
