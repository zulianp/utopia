#ifndef UTOPIA_MOONOLITH_STK_HPP
#define UTOPIA_MOONOLITH_STK_HPP

#include "utopia_ConvertMesh.hpp"
#include "utopia_moonolith_ForwardDeclarations.hpp"
#include "utopia_stk_ForwardDeclarations.hpp"

namespace utopia {
    template <>
    class ConvertMesh<utopia::stk::Mesh, utopia::moonolith::Mesh> {
    public:
        static void apply(const utopia::stk::Mesh &in, utopia::moonolith::Mesh &out);
    };

    // Eventually implement this
    // template <>
    // class ConvertMesh<utopia::moonolith::Mesh, utopia::stk::Mesh> {
    // public:
    //     static void apply(const utopia::moonolith::Mesh &in, utopia::stk::Mesh &out);
    // };
}  // namespace utopia

#endif  // UTOPIA_MOONOLITH_STK_HPP
