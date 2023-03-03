

#include "utopia_Base.hpp"
#include "utopia_polymorphic_QPSolver_impl.hpp"

#ifdef UTOPIA_ENABLE_PETSC
#include "utopia_petsc_Types.hpp"
#endif  // UTOPIA_ENABLE_PETSC

namespace utopia {

#ifdef UTOPIA_ENABLE_PETSC
    template class OmniQPSolver<PetscMatrix, PetscVector>;
#endif  // UTOPIA_ENABLE_PETSC

}  // namespace utopia