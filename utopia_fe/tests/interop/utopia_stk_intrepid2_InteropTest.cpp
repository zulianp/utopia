#include "utopia_fe_base.hpp"

#ifdef UTOPIA_WITH_INTREPID2

#include "utopia_Testing.hpp"
#include "utopia_ui.hpp"

#include "utopia_stk_FunctionSpace.hpp"
#include "utopia_stk_Mesh.hpp"

#include "utopia_FEInteroperability.hpp"

void poisson_problem() {}

void interop_stk_intrepid2() { UTOPIA_RUN_TEST(poisson_problem); }

UTOPIA_REGISTER_TEST_FUNCTION(interop_stk_intrepid2);

#endif  // UTOPIA_WITH_INTREPID2