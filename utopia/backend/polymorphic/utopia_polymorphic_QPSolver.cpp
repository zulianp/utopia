

#include "utopia_Base.hpp"
#include "utopia_polymorphic_QPSolver_impl.hpp"

#ifdef UTOPIA_WITH_PETSC
#include "utopia_petsc_Types.hpp"
#endif  // UTOPIA_WITH_PETSC

namespace utopia {

#ifdef UTOPIA_WITH_PETSC
    template class PolymorphicQPSolver<PetscMatrix, PetscVector>;
#endif  // UTOPIA_WITH_PETSC

}  // namespace utopia