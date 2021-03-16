#include "utopia_Testing.hpp"

#include "utopia_moonolith_Mesh.hpp"
#include "utopia_stk_Mesh.hpp"
#include "utopia_ui.hpp"

#include "utopia_FEInteroperability.hpp"

using namespace utopia;

void stk_moonolith_convert() {
    using MeshFrom = utopia::stk::Mesh;
    using MeshTo = utopia::moonolith::Mesh;

    MeshFrom mesh_from;
    utopia_test_assert(mesh_from.read("../data/knf/rectangle_4_tris.e"));

    MeshTo mesh_to;
    convert_mesh(mesh_from, mesh_to);
    utopia_test_assert(mesh_to.write("membrane.vtu"));
}

void interop() { UTOPIA_RUN_TEST(stk_moonolith_convert); }

UTOPIA_REGISTER_TEST_FUNCTION(interop);
