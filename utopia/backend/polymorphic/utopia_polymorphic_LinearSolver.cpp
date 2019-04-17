#include "utopia_polymorphic_LinearSolver_impl.hpp"
#include "utopia_Base.hpp"

#ifdef WITH_PETSC
#include "utopia_petsc_Types.hpp"
#include "utopia_petsc_LinearSolverFactory.hpp"
#endif

namespace utopia {
#ifdef WITH_PETSC
    template class PolymorphicLinearSolver<DSMatrixd, DVectord>;
#endif
}
