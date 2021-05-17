#ifndef UTOPIA_MOONOLITH_STK_HPP
#define UTOPIA_MOONOLITH_STK_HPP

#include "utopia_ConvertFunctionSpace.hpp"
#include "utopia_ConvertMesh.hpp"
#include "utopia_ExtractSurface.hpp"
#include "utopia_ExtractTraceSpace.hpp"

#include "utopia_moonolith_ForwardDeclarations.hpp"
#include "utopia_stk_ForwardDeclarations.hpp"

namespace utopia {
    template <>
    class ConvertMesh<utopia::stk::Mesh, utopia::moonolith::Mesh> {
    public:
        static void apply(const utopia::stk::Mesh &in, utopia::moonolith::Mesh &out);
    };

    template <>
    class ConvertFunctionSpace<utopia::stk::FunctionSpace, utopia::moonolith::FunctionSpace> {
    public:
        static void apply(const utopia::stk::FunctionSpace &in, utopia::moonolith::FunctionSpace &out);
    };

    template <>
    class ExtractSurface<utopia::stk::Mesh, utopia::moonolith::Mesh> {
    public:
        static void apply(const utopia::stk::Mesh &in, utopia::moonolith::Mesh &out);
    };

    template <>
    class ExtractTraceSpace<utopia::stk::FunctionSpace, utopia::moonolith::FunctionSpace> {
    public:
        static void apply(const utopia::stk::FunctionSpace &in, utopia::moonolith::FunctionSpace &out);
    };

}  // namespace utopia

#endif  // UTOPIA_MOONOLITH_STK_HPP
