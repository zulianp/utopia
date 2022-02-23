#include "utopia_Base.hpp"
#include "utopia_polymorphic_LinearSolver_impl.hpp"

#ifdef UTOPIA_WITH_PETSC
#include "utopia_petsc_LinearSolverFactory.hpp"
#include "utopia_petsc_Types.hpp"
#endif

#ifdef UTOPIA_WITH_TRILINOS
#include "utopia_trilinos_LinearSolverFactory.hpp"
#include "utopia_trilinos_Types.hpp"
#endif

namespace utopia {
#ifdef UTOPIA_WITH_PETSC
    template class OmniLinearSolver<PetscMatrix, PetscVector>;
#endif
#ifdef UTOPIA_WITH_TRILINOS
    template class OmniLinearSolver<TpetraMatrix, TpetraVector>;
#endif
}  // namespace utopia
