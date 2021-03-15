#include "utopia_Testing.hpp"

#include "utopia_moonolith_Mesh.hpp"
#include "utopia_ui.hpp"

using namespace utopia;

void moonolith_create_mesh() {
    using Mesh_t = utopia::moonolith::Mesh;

    Mesh_t mesh;
    utopia_test_assert(mesh.read("/Users/zulianp/Desktop/code/par_moonolith/examples/data/triangle_105.tri"));
    utopia_test_assert(mesh.write("prova.vtu"));
}

void umoonolith() { UTOPIA_RUN_TEST(moonolith_create_mesh); }

UTOPIA_REGISTER_TEST_FUNCTION(umoonolith);
