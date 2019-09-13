#include "utopia_petsc_SNES_impl.hpp"
#include "utopia_petsc_Types.hpp"

namespace utopia {
    template class SNESSolver<DSMatrixd, DVectord, PETSC>;
    //FIXME
    // template class SNESSolver<DMatrixd, DVectord, PETSC>;
}
