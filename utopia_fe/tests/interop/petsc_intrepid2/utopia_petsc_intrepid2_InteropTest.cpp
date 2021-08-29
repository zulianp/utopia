#include "utopia_fe_base.hpp"

#ifdef UTOPIA_WITH_PETSC

#include "utopia_Testing.hpp"
#include "utopia_ui.hpp"

#include "utopia_FEInteroperability.hpp"

// Stk includes
#include "utopia_petsc_FunctionSpace.hpp"
#include "utopia_petsc_Mesh.hpp"

// utopia/kokkos includes
#include "utopia_kokkos_L2Norm.hpp"
#include "utopia_kokkos_LaplaceOperator.hpp"
#include "utopia_kokkos_LinearElasticity.hpp"
#include "utopia_kokkos_NeoHookean.hpp"
#include "utopia_kokkos_VectorLaplaceOperator.hpp"

// utopia/intrepid2 includes
#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_ShellTools.hpp"

// Interop includes
#include "utopia_petsc_intrepid2.hpp"

// FIXME: This is the last include because the operator files are not yet in the correct place
#include "utopia_UnitCubeSpaceAndFETest.hpp"

// #include "utopia_ScalarProductTest.hpp"

using namespace utopia;

using PetscScalar_t = Traits<utopia::petsc::FunctionSpace>::Scalar;
void interop_petsc_intrepid2() {
    UnitCubeSpaceAndFETest<utopia::petsc::FunctionSpace, utopia::intrepid2::FE<PetscScalar_t>>().run();
}

UTOPIA_REGISTER_TEST_FUNCTION(interop_petsc_intrepid2);

// void scalar_product_petsc_intrepid2() { ScalarProductTest<utopia::petsc::FunctionSpace>().run(); }

// UTOPIA_REGISTER_TEST_FUNCTION(scalar_product_petsc_intrepid2);

#endif  // UTOPIA_WITH_PETSC
