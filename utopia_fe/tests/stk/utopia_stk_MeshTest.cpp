#include "utopia_Testing.hpp"

#include "utopia_stk_Mesh.hpp"
#include "utopia_ui.hpp"

using namespace utopia;

void stk_create_mesh() {
    using Mesh_t = utopia::stk::Mesh;

    InputParameters params;
    params.set("type", "file");
    params.set("path", "../data/knf/pump/membrane.e");

    Mesh_t mesh;
    mesh.read(params);
}

void ustk() { UTOPIA_RUN_TEST(stk_create_mesh); }

UTOPIA_REGISTER_TEST_FUNCTION(ustk);
