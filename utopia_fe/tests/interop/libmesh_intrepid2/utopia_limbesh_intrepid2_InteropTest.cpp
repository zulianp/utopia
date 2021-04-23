#include "utopia_fe_base.hpp"

#ifdef UTOPIA_WITH_INTREPID2

#include "utopia_Testing.hpp"
#include "utopia_ui.hpp"

#include "utopia_FEInteroperability.hpp"

// Stk includes
#include "utopia_libmesh_FunctionSpace_new.hpp"
#include "utopia_libmesh_Mesh.hpp"

// Intrepid2 includes
#include "utopia_intrepid2_LaplaceOperator.hpp"
#include "utopia_intrepid2_LinearElasticity.hpp"
#include "utopia_intrepid2_ShellTools.hpp"
#include "utopia_intrepid2_VectorLaplaceOperator.hpp"

#include "utopia_libmesh_intrepid2.hpp"

// FIXME: This is the last include because the operator files are not yet in the correct place
#include "utopia_SpaceAndFETest.hpp"

using namespace utopia;

using StkScalar = Traits<utopia::libmesh::FunctionSpace>::Scalar;
void interop_libmesh_intrepid2() {
    SpaceAndFETest<utopia::libmesh::FunctionSpace, utopia::intrepid2::FE<StkScalar>>().run();
}

UTOPIA_REGISTER_TEST_FUNCTION(interop_libmesh_intrepid2);

#endif  // UTOPIA_WITH_INTREPID2
