#include "utopia_petsc_TaoQPSolver_impl.hpp"

#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_Vector_impl.hpp"

namespace utopia {
    template class TaoQPSolver<PetscMatrix, PetscVector>;
}  // namespace utopia
