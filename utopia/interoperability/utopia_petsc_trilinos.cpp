#include "utopia_petsc_trilinos.hpp"

#ifdef WITH_TRILINOS
#ifdef WITH_PETSC

#include "utopia_trilinos.hpp"

namespace utopia {
	template class KSPSolver<TSMatrixd, TVectord, TRILINOS>;
	template class Factorization<TSMatrixd, TVectord, TRILINOS>;
}

#endif //WITH_PETSC
#endif //WITH_TRILINOS
