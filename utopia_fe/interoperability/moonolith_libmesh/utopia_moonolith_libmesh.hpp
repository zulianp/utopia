#ifndef UTOPIA_MOONOLITH_LIBMESH_HPP
#define UTOPIA_MOONOLITH_LIBMESH_HPP

// Templates
#include "utopia_ConvertFunctionSpace.hpp"
#include "utopia_ConvertMesh.hpp"
#include "utopia_ExtractTraceSpace.hpp"

// Forward type
#include "utopia_libmesh_ForwardDeclarations.hpp"
#include "utopia_moonolith_ForwardDeclarations.hpp"

namespace utopia {
    template <>
    class ConvertMesh<utopia::libmesh::Mesh, utopia::moonolith::Mesh> {
    public:
        static void apply(const utopia::libmesh::Mesh &in, utopia::moonolith::Mesh &out);
    };

    template <>
    class ConvertFunctionSpace<utopia::libmesh::FunctionSpace, utopia::moonolith::FunctionSpace> {
    public:
        static void apply(const utopia::libmesh::FunctionSpace &in, utopia::moonolith::FunctionSpace &out);
    };

    template <>
    class ExtractTraceSpace<utopia::libmesh::FunctionSpace, utopia::moonolith::FunctionSpace> {
    public:
        static void apply(const utopia::libmesh::FunctionSpace &in, utopia::moonolith::FunctionSpace &out);
    };

}  // namespace utopia

#endif  // UTOPIA_MOONOLITH_LIBMESH_HPP