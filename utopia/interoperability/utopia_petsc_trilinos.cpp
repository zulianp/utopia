#include "utopia_petsc_trilinos.hpp"
#include "utopia_trilinos.hpp"

namespace utopia {
	template class KSPSolver<TSMatrixd, TVectord, TRILINOS>;
}