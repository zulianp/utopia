#include "utopia_petsc_DMPlex_FunctionSpace.hpp"
#include "utopia_Quad4.hpp"

namespace utopia {
    template class FunctionSpace<PetscDMPlex<StaticVector<PetscScalar, 2>, ArrayView<PetscInt, 2>>, 1>;
}
