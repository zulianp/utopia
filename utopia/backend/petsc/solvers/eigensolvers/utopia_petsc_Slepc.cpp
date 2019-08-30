#include "utopia_petsc_Slepc_impl.hpp"
#include "utopia_petsc_Types.hpp"

namespace utopia {
    template class SlepcSolver<DSMatrixd, DVectord, PETSC_EXPERIMENTAL>;
    template class SlepcSolver<DMatrixd, DVectord, PETSC_EXPERIMENTAL>;
}
