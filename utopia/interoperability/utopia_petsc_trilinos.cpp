#include "utopia_petsc_trilinos.hpp"

#ifdef WITH_TRILINOS
#ifdef WITH_PETSC

#include "utopia_trilinos.hpp"

namespace utopia {
    template class KSPSolver<TpetraMatrixd, TpetraVectord, TRILINOS>;
    template class Factorization<TpetraMatrixd, TpetraVectord, TRILINOS>;
    template class BiCGStab<TpetraMatrixd, TpetraVectord, TRILINOS>;
    template class GaussSeidel<TpetraMatrixd, TpetraVectord, TRILINOS>;
}  // namespace utopia

#endif //WITH_PETSC
#endif //WITH_TRILINOS
