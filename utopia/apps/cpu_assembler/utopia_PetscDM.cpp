#include "utopia_PetscDM_impl.hpp"

namespace utopia {
    template class PetscNode<1>;
    template class PetscNode<2>;
    template class PetscNode<3>;

    template class PetscDM<1>;
    template class PetscDM<2>;
    template class PetscDM<3>;

    template class FunctionSpace<PetscDM<2>, 1, PetscUniformQuad4>;
    template class FunctionSpace<PetscDM<3>, 1, PetscUniformHex8>;
}
