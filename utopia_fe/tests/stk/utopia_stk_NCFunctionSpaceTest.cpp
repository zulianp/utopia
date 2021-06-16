
#include "utopia_Testing.hpp"

#include "utopia_NCFunctionSpaceTest.hpp"

#include "utopia_moonolith_stk_FETransfer.hpp"
#include "utopia_stk_FunctionSpace.hpp"
#include "utopia_stk_Mesh.hpp"

#include "utopia_ui.hpp"

#include "utopia_RunParallelTest.hpp"

using namespace utopia;

void stk_ncfunctionspace() { utopia::run_parallel_test<NCFunctionSpaceTest<utopia::stk::FunctionSpace>>(); }

UTOPIA_REGISTER_TEST_FUNCTION(stk_ncfunctionspace);
