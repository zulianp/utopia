#include "utopia_petsc_SemismoothNewton.hpp"
#include "utopia_petsc_SemismoothNewton_impl.hpp"

#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_Vector_impl.hpp"

namespace utopia {
    template class SemismoothNewton<PetscMatrix, PetscVector, PETSC_EXPERIMENTAL>;
}  // namespace utopia
