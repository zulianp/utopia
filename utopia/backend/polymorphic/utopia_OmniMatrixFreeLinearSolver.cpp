#include "utopia_Base.hpp"
#include "utopia_OmniMatrixFreeLinearSolver_impl.hpp"

#ifdef UTOPIA_ENABLE_PETSC
#include "utopia_petsc_LinearSolverFactory.hpp"
#include "utopia_petsc_Types.hpp"
#endif

#ifdef UTOPIA_ENABLE_TRILINOS
#include "utopia_trilinos_LinearSolverFactory.hpp"
#include "utopia_trilinos_Types.hpp"
#endif

namespace utopia {
#ifdef UTOPIA_ENABLE_PETSC
    template class OmniMatrixFreeLinearSolver<PetscVector>;
#endif
#ifdef UTOPIA_ENABLE_TRILINOS
    template class OmniMatrixFreeLinearSolver<TpetraVector>;
#endif
}  // namespace utopia
