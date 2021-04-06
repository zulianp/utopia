
#include "utopia_Testing.hpp"

#include "utopia_MeshTest.hpp"

#include "utopia_stk_Mesh.hpp"
#include "utopia_ui.hpp"

#include "utopia_RunParallelTest.hpp"

using namespace utopia;

void ustk() { utopia::run_parallel_test<MeshTest<utopia::stk::Mesh>>(); }

UTOPIA_REGISTER_TEST_FUNCTION(ustk);
