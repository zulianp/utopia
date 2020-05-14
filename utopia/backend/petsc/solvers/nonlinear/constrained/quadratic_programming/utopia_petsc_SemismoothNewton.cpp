#include "utopia_petsc_SemismoothNewton.hpp"
#include "utopia_petsc_SemismoothNewton_impl.hpp"

namespace utopia {
    template class SemismoothNewton<PetscMatrix, PetscVector, PETSC_EXPERIMENTAL>;
}  // namespace utopia
