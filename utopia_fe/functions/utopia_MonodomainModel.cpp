#include "utopia_MonodomainModel.hpp"
#include "utopia_fe_base.hpp"

namespace utopia {
	template class MonodomainModel<LibMeshFunctionSpace, USMatrix, UVector>;
}