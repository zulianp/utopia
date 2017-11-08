#ifndef UTOPIA_LIBMESH_FUNCTIONAL_TRAITS_HPP
#define UTOPIA_LIBMESH_FUNCTIONAL_TRAITS_HPP 

#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_libmesh_AssemblyContext.hpp"

namespace utopia {
	template<>
	class FunctionalTraits<LibMeshFunctionSpace, AssemblyContext<LIBMESH_TAG> > {
	public:
		inline static int type(const LibMeshFunctionSpace &space, const AssemblyContext<LIBMESH_TAG> &ctx)
		{
			return utopia::POLYNOMIAL_FUNCTION;
		}

		inline static int order(const LibMeshFunctionSpace &space, const AssemblyContext<LIBMESH_TAG> &ctx)
		{
			return space.order(ctx.current_element());
		}
	};
}

#endif //UTOPIA_LIBMESH_FUNCTIONAL_TRAITS_HPP
