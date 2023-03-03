#include "utopia_petsc_trilinos.hpp"

#ifdef UTOPIA_WITH_TRILINOS
#ifdef UTOPIA_ENABLE_PETSC

#include <Kokkos_Macros.hpp>
#ifndef KOKKOS_ENABLE_CUDA
#include "utopia_trilinos.hpp"

namespace utopia {
    template class KSPSolver<TpetraMatrixd, TpetraVectord, TRILINOS>;
    template class Factorization<TpetraMatrixd, TpetraVectord, TRILINOS>;
    template class BiCGStab<TpetraMatrixd, TpetraVectord, TRILINOS>;
    template class GaussSeidel<TpetraMatrixd, TpetraVectord, TRILINOS>;
}  // namespace utopia

#endif  // not def KOKKOS_ENABLE_CUDA
#endif  // UTOPIA_ENABLE_PETSC
#endif  // UTOPIA_WITH_TRILINOS
