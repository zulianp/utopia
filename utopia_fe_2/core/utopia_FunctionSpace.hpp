#ifndef UTOPIA_FUNCTION_SPACE_HPP
#define UTOPIA_FUNCTION_SPACE_HPP

#include "utopia_Mesh.hpp"

namespace utopia {

    class IFunctionSpace : public Configurable, public Describable {};

    template <class... T>
    class FunctionSpace {};
}  // namespace utopia

#endif  // UTOPIA_FUNCTION_SPACE_HPP
