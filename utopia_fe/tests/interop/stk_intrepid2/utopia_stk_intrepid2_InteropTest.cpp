#include "utopia_fe_base.hpp"

#ifdef UTOPIA_WITH_INTREPID2

#include "utopia_Testing.hpp"
#include "utopia_ui.hpp"

#include "utopia_stk_FunctionSpace.hpp"
#include "utopia_stk_Mesh.hpp"

#include "utopia_FEInteroperability.hpp"

#include "utopia_intrepid2_LaplaceOperator.hpp"
#include "utopia_intrepid2_LinearElasticity.hpp"
// #include "utopia_intrepid2_ShellFE.hpp"
#include "utopia_intrepid2_ShellTools.hpp"
#include "utopia_intrepid2_VectorLaplaceOperator.hpp"

#include "utopia_SpaceAndFETest.hpp"

using namespace utopia;

using StkScalar = Traits<utopia::stk::FunctionSpace>::Scalar;
void interop_stk_intrepid2() { SpaceAndFETest<utopia::stk::FunctionSpace, utopia::intrepid2::FE<StkScalar>>().run(); }

UTOPIA_REGISTER_TEST_FUNCTION(interop_stk_intrepid2);

#endif  // UTOPIA_WITH_INTREPID2
