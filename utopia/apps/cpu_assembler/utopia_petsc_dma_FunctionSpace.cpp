#include "utopia_petsc_dma_FunctionSpace_impl.hpp"
#include "utopia_petsc_DirichletBoundaryConditions.hpp"
#include "utopia_Tri3.hpp"

namespace utopia {
    template class FunctionSpace<PetscDM<2>, 1, PetscUniformQuad4>;
    template class FunctionSpace<PetscDM<3>, 1, PetscUniformHex8>;

    template class FunctionSpace<PetscDM<2>, 2, PetscUniformQuad4>;
    template class FunctionSpace<PetscDM<3>, 3, PetscUniformHex8>;

    template class FunctionSpace<PetscDM<2>, 3, PetscUniformQuad4>;
    template class FunctionSpace<PetscDM<3>, 4, PetscUniformHex8>;


    /////////////////////////////////////////////////////////////////

    template class FunctionSpace<PetscDM<2>, 1, Tri3<PetscScalar, 2>>;
    template class FunctionSpace<PetscDM<2>, 2, Tri3<PetscScalar, 2>>;
}
