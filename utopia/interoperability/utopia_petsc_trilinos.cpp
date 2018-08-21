#include "utopia_petsc_trilinos.hpp"

#ifdef WITH_TRILINOS

#include "utopia_trilinos.hpp"

namespace utopia {
	template class KSPSolver<TSMatrixd, TVectord, TRILINOS>;
}

#endif //WITH_TRILINOS
