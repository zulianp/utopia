#include "utopia_Testing.hpp"

#include "utopia_stk_Mesh.hpp"
#include "utopia_ui.hpp"

using namespace utopia;

void stk_create_mesh() {
    using Mesh_t = utopia::stk::Mesh;

    Mesh_t mesh;
    utopia_test_assert(mesh.read("../data/knf/pump/membrane.e"));
    utopia_test_assert(mesh.write("prova.e"));
}

void ustk() { UTOPIA_RUN_TEST(stk_create_mesh); }

UTOPIA_REGISTER_TEST_FUNCTION(ustk);
