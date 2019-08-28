#include "utopia_petsc_Eval_Cond_impl.hpp"

#ifdef WITH_SLEPC
#include "utopia_petsc.hpp"

namespace utopia {
    template class Cond<DSMatrixd, PETSC>;
    template class Cond<DMatrixd, PETSC>;
}

#endif //WITH_SLEPC
