#include "utopia_petsc_Eval_Cond_impl.hpp"

#ifdef UTOPIA_ENABLE_SLEPC
#include "utopia_petsc_Matrix.hpp"

namespace utopia {
    // template class Cond<PetscMatrix, PETSC>;
    // template class Cond<PetscMatrix, PETSC>;

    template class Cond<PetscMatrix, PETSC>;
}  // namespace utopia

#endif  // UTOPIA_ENABLE_SLEPC
