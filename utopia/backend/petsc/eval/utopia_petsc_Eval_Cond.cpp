#include "utopia_petsc_Eval_Cond_impl.hpp"

#ifdef WITH_SLEPC
#include "utopia_petsc_Matrix.hpp"

namespace utopia {
    // template class Cond<PetscMatrix, PETSC>;
    // template class Cond<PetscMatrix, PETSC>;

    template class Cond<PetscMatrix, PETSC>;
}  // namespace utopia

#endif  // WITH_SLEPC
