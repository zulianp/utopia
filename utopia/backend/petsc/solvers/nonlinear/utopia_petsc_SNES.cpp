#include "utopia_petsc_SNES_impl.hpp"
#include "utopia_petsc_Types.hpp"

#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_Vector_impl.hpp"

namespace utopia {
    template class SNESSolver<PetscMatrix, PetscVector, PETSC>;
    // FIXME
    // template class SNESSolver<PetscMatrix, PetscVector, PETSC>;
}  // namespace utopia
