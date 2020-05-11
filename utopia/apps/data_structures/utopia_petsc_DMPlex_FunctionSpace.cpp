#include "utopia_petsc_Base.hpp"

#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 12, 4)

#include "utopia_Quad4.hpp"
#include "utopia_Tri3.hpp"
#include "utopia_petsc_DMPlex_FunctionSpace.hpp"

namespace utopia {
    // template class FunctionSpace<PetscDMPlex<StaticVector<PetscScalar, 2>, ArrayView<PetscInt, 4>>,
    //                              1,
    //                              Quad4<PetscScalar, 2>>;

    template class FunctionSpace<PetscDMPlex<StaticVector<PetscScalar, 2>, ArrayView<PetscInt, 3>>,
                                 1,
                                 Tri3<PetscScalar, 2>>;

}  // namespace utopia

#endif  // UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 12, 4)
