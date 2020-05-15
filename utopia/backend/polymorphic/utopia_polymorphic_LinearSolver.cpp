#include "utopia_Base.hpp"
#include "utopia_polymorphic_LinearSolver_impl.hpp"

#ifdef WITH_PETSC
#include "utopia_petsc_LinearSolverFactory.hpp"
#include "utopia_petsc_Types.hpp"
#endif

namespace utopia {
#ifdef WITH_PETSC
    template class PolymorphicLinearSolver<PetscMatrix, PetscVector>;
#endif
}  // namespace utopia
