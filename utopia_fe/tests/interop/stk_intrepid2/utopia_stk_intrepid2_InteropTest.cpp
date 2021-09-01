#include "utopia_fe_base.hpp"

#ifdef UTOPIA_WITH_INTREPID2

#include "utopia_Testing.hpp"
#include "utopia_ui.hpp"

#include "utopia_FEInteroperability.hpp"

// Stk includes
#include "utopia_stk_FunctionSpace.hpp"
#include "utopia_stk_Mesh.hpp"

// utopia/kokkos includes
#include "utopia_kokkos_L2Norm.hpp"
#include "utopia_kokkos_LaplaceOperator.hpp"
#include "utopia_kokkos_LinearElasticity.hpp"
#include "utopia_kokkos_NeoHookean.hpp"
#include "utopia_kokkos_VectorLaplaceOperator.hpp"

// utopia/intrepid2 includes
#include "utopia_intrepid2_ShellTools.hpp"

// Interop includes
#include "utopia_stk_intrepid2.hpp"

// FIXME: This is the last include because the operator files are not yet in the correct place
#include "utopia_SpaceAndFETest.hpp"

#include "utopia_ScalarProductTest.hpp"

using namespace utopia;

using StkScalar = Traits<utopia::stk::FunctionSpace>::Scalar;
void interop_stk_intrepid2() { SpaceAndFETest<utopia::stk::FunctionSpace, utopia::intrepid2::FE<StkScalar>>().run(); }

UTOPIA_REGISTER_TEST_FUNCTION(interop_stk_intrepid2);

void scalar_product_stk_intrepid2() { ScalarProductTest<utopia::stk::FunctionSpace>().run(); }

UTOPIA_REGISTER_TEST_FUNCTION(scalar_product_stk_intrepid2);

#endif  // UTOPIA_WITH_INTREPID2
