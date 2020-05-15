#include "utopia_petsc_SNES_impl.hpp"
#include "utopia_petsc_Types.hpp"

namespace utopia {
    template class SNESSolver<PetscMatrix, PetscVector, PETSC>;
    // FIXME
    // template class SNESSolver<PetscMatrix, PetscVector, PETSC>;
}  // namespace utopia
