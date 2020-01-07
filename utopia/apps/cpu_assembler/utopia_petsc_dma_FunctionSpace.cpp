#include "utopia_petsc_dma_FunctionSpace_impl.hpp"

namespace utopia {
    template class FunctionSpace<PetscDM<2>, 1, PetscUniformQuad4>;
    template class FunctionSpace<PetscDM<3>, 1, PetscUniformHex8>;
}
