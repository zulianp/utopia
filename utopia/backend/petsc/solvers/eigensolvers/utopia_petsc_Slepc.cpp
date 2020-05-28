#include "utopia_petsc_Slepc_impl.hpp"
#include "utopia_petsc_Types.hpp"

#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_Vector_impl.hpp"

namespace utopia {
    // template class SlepcSolver<PetscMatrix, PetscVector, PETSC_EXPERIMENTAL>;
    // template class SlepcSolver<PetscMatrix, PetscVector, PETSC_EXPERIMENTAL>;

    template class SlepcSolver<PetscMatrix, PetscVector, PETSC_EXPERIMENTAL>;
}  // namespace utopia
