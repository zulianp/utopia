#include "utopia_m3elinsol_impl.hpp"

#ifdef UTOPIA_WITH_TRILINOS
#include "utopia_trilinos.hpp"
#endif  // UTOPIA_WITH_TRILINOS

#ifdef UTOPIA_WITH_PETSC
#include "utopia_petsc.hpp"
#include "utopia_petsc_impl.hpp"
#endif  // UTOPIA_WITH_PETSC

#ifdef UTOPIA_WITH_BLAS
#include "utopia_blas.hpp"
#endif  // UTOPIA_WITH_BLAS

// explicit instantiations

namespace utopia {

#ifdef UTOPIA_WITH_PETSC
    template class ASPAMG<PetscMatrix, PetscVector>;
#endif  // UTOPIA_WITH_PETSC

#ifdef UTOPIA_WITH_TRILINOS

    template class ASPAMG<TpetraMatrixd, TpetraVectord>;
#endif  // UTOPIA_WITH_TRILINOS

}  // namespace utopia
