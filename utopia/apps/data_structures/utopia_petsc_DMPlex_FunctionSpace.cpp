#include "utopia_petsc_Base.hpp"

#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 12, 4)

#include "utopia_Quad4.hpp"
#include "utopia_petsc_DMPlex_FunctionSpace.hpp"

namespace utopia {
    template class FunctionSpace<PetscDMPlex<StaticVector<PetscScalar, 2>, ArrayView<PetscInt, 2>>, 1>;
}

#endif  // UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 12, 4)
