

#include "utopia_polymorphic_QPSolver_impl.hpp"
#include "utopia_Base.hpp"

#ifdef WITH_PETSC
#include "utopia_petsc_Types.hpp"
#endif //WITH_PETSC

namespace utopia {

#ifdef WITH_PETSC
        template class PolymorphicQPSolver<PetscMatrix, PetscVector>;
#endif //WITH_PETSC

}