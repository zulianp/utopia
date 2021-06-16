
#include "utopia_Testing.hpp"

#include "utopia_NCFunctionSpaceTest.hpp"

#include "utopia_libmesh_FunctionSpace_new.hpp"
#include "utopia_libmesh_Mesh.hpp"
#include "utopia_moonolith_libmesh_FETransfer.hpp"

#include "utopia_ui.hpp"

#include "utopia_RunParallelTest.hpp"

using namespace utopia;

void libmesh_ncfunctionspace() { utopia::run_parallel_test<NCFunctionSpaceTest<utopia::libmesh::FunctionSpace>>(); }

UTOPIA_REGISTER_TEST_FUNCTION(libmesh_ncfunctionspace);
