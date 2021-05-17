#ifndef UTOPIA_LIBMESH_FUNCTIONAL_TRAITS_HPP
#define UTOPIA_LIBMESH_FUNCTIONAL_TRAITS_HPP

#include "utopia_FunctionalTraits.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"

// #include "utopia_libmesh_AssemblyContext.hpp"

namespace utopia {
    template <class AssemblyContextT>
    class FunctionalTraits<LibMeshFunctionSpace, AssemblyContextT> {
    public:
        inline static int type(const LibMeshFunctionSpace &space, const AssemblyContextT &ctx) {
            return utopia::POLYNOMIAL_FUNCTION;
        }

        inline static int order(const LibMeshFunctionSpace &space, const AssemblyContextT &ctx) {
            return space.order(ctx.current_element());
        }
    };
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_FUNCTIONAL_TRAITS_HPP
