#ifndef UTOPIA_SFEM_STK_EXTRACTSURFACE_HPP
#define UTOPIA_SFEM_STK_EXTRACTSURFACE_HPP

#include "utopia_ExtractSurface.hpp"

#include "utopia_sfem_ForwardDeclarations.hpp"
#include "utopia_stk_ForwardDeclarations.hpp"

namespace utopia {
    template <>
    class ExtractSurface<utopia::stk::Mesh, utopia::sfem::Mesh> {
    public:
        static void apply(const utopia::stk::Mesh &in, utopia::sfem::Mesh &out);
    };
}  // namespace utopia

#endif  // UTOPIA_SFEM_STK_EXTRACTSURFACE_HPP