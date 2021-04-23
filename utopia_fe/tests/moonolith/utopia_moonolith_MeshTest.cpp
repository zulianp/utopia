#include "utopia_Testing.hpp"

#include "utopia_moonolith_FETransfer.hpp"
#include "utopia_moonolith_FunctionSpace.hpp"
#include "utopia_moonolith_Mesh.hpp"
#include "utopia_ui.hpp"

using namespace utopia;

static Path db_path = "../data/testing";

void moonolith_create_mesh() {
    using Mesh_t = utopia::moonolith::Mesh;

    Mesh_t mesh;
    utopia_test_assert(mesh.read(db_path / "triangle_105.tri"));
    utopia_test_assert(mesh.write("prova.vtu"));
}

void moonolith_create_space() {
    using Space_t = utopia::moonolith::FunctionSpace;

    Space_t space;
    auto params = param_list(param("path", db_path / "triangle_105.tri"));
    space.read(params);

    utopia_test_assert(space.mesh().write("prova.vtu"));
}

void moonolith_fe_transfer() {
    using Space_t = utopia::moonolith::FunctionSpace;

    Space_t space_from;
    auto params_from = param_list(param("path", db_path / "triangle_4.tri"));
    space_from.read(params_from);

    auto params_to = param_list(param("path", db_path / "square_200.tri"));
    Space_t space_to;
    space_to.read(params_to);

    utopia::moonolith::FETransfer transfer;
    utopia_test_assert(transfer.init(make_ref(space_from), make_ref(space_to)));
}

void umoonolith() {
    if (mpi_world_size() == 1) {
        // Only serial
        UTOPIA_RUN_TEST(moonolith_create_mesh);
        UTOPIA_RUN_TEST(moonolith_create_space);
        UTOPIA_RUN_TEST(moonolith_fe_transfer);
    }
}

UTOPIA_REGISTER_TEST_FUNCTION(umoonolith);
