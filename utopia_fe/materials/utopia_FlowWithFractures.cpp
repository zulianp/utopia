#include "utopia_FlowWithFractures.hpp"
#include "utopia_libmesh.hpp"

namespace utopia {
    template class FlowWithFractures<LibMeshFunctionSpace, USparseMatrix, UVector>;
}
