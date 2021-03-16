#include "utopia_Testing.hpp"

#include "utopia_moonolith_FunctionSpace.hpp"
#include "utopia_moonolith_Mesh.hpp"
#include "utopia_ui.hpp"

using namespace utopia;

void moonolith_create_mesh() {
    using Mesh_t = utopia::moonolith::Mesh;

    Mesh_t mesh;
    utopia_test_assert(mesh.read("/Users/zulianp/Desktop/code/par_moonolith/examples/data/triangle_105.tri"));
    utopia_test_assert(mesh.write("prova.vtu"));
}

void moonolith_create_space() {
    using Space_t = utopia::moonolith::FunctionSpace;

    Space_t space;
    auto params = param_list(param("path", "/Users/zulianp/Desktop/code/par_moonolith/examples/data/triangle_105.tri"));
    space.read(params);

    utopia_test_assert(space.mesh().write("prova.vtu"));
}

void umoonolith() {
    UTOPIA_RUN_TEST(moonolith_create_mesh);
    UTOPIA_RUN_TEST(moonolith_create_space);
}

UTOPIA_REGISTER_TEST_FUNCTION(umoonolith);
