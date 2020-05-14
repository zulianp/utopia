#include "utopia_petsc_DMDA_FunctionSpace.hpp"
#include "utopia_UniformQuad4.hpp"

namespace utopia {
    template class FunctionSpace< PetscDMDA<StaticVector<PetscScalar, 2>, ArrayView<PetscInt, 2>>, 1, UniformQuad4<PetscScalar>>;
}  // namespace utopia
