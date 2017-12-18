#ifndef UTOPIA_LIBMESH_FE_FORWARD_DECLARATIONS_HPP
#define UTOPIA_LIBMESH_FE_FORWARD_DECLARATIONS_HPP 

#include "utopia_FEForwardDeclarations.hpp"

#ifdef NDEBUG
#if defined(DEBUG)
#error "DEBUG macro defined when it should not be defined"
#endif
#endif //NDEBUG

static const int LIBMESH_TAG = 1001;

namespace utopia {
	class LibMeshFunctionSpace;
}

#endif //UTOPIA_LIBMESH_FE_FORWARD_DECLARATIONS_HPP
