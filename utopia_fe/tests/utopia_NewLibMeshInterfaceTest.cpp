#include "utopia_Testing.hpp"

#include "utopia_ui.hpp"

#include "utopia_libmesh_Mesh.hpp"

using namespace utopia;

void lm_create_mesh() {
    using Mesh_t = utopia::libmesh::Mesh;

    InputParameters params;
    params.set("type", "cube");
    params.set("elem_type", "TET4");
    params.set("show", true);

    Mesh_t mesh;
    mesh.read(params);
}

void lm() { UTOPIA_RUN_TEST(lm_create_mesh); }

UTOPIA_REGISTER_TEST_FUNCTION(lm);
