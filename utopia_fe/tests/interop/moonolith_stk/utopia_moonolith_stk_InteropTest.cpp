#include "utopia_Testing.hpp"

#include "utopia_moonolith_Mesh.hpp"
#include "utopia_stk_Mesh.hpp"

#include "utopia_moonolith_FunctionSpace.hpp"
#include "utopia_moonolith_stk_FETransfer.hpp"
#include "utopia_stk_FunctionSpace.hpp"

#include "utopia_ui.hpp"

#include "utopia_FEInteroperability.hpp"

using namespace utopia;

static const int nx = 10, ny = 10, nz = 10;

void stk_moonolith_convert_mesh() {
    using MeshFrom = utopia::stk::Mesh;
    using MeshTo = utopia::moonolith::Mesh;

    MeshFrom mesh_from;
    mesh_from.unit_cube(nx, ny, nz);

    MeshTo mesh_to;
    convert_mesh(mesh_from, mesh_to);
    // utopia_test_assert(mesh_to.write("dump.vtu"));
}

void stk_moonolith_extract_surface() {
    using MeshFrom = utopia::stk::Mesh;
    using MeshTo = utopia::moonolith::Mesh;

    MeshFrom volume;
    volume.unit_cube(nx, ny, nz);

    MeshTo surface;
    extract_surface(volume, surface);

    if (volume.comm().size() == 1) {
        // surface.write("surf.vtu");
        utopia_test_assert(surface.n_nodes() == (ny + 1) * (nz + 1) * 2);
        utopia_test_assert(surface.manifold_dimension() == 2);
        utopia_test_assert(surface.spatial_dimension() == 3);
    }
}

void stk_moonolith_convert_space() {
    using FunctionSpaceFrom = utopia::stk::FunctionSpace;
    using FunctionSpaceTo = utopia::moonolith::FunctionSpace;

    auto params = param_list(param("mesh", param_list(param("type", "cube"))));
    FunctionSpaceFrom space_from;
    space_from.read(params);

    FunctionSpaceTo space_to;
    convert_function_space(space_from, space_to);

    // disp(space_to.n_dofs());
    utopia_test_assert(space_to.mesh().write("dump.vtu"));
}

void stk_moonolith_extract_trace_space() {
    using FunctionSpaceFrom = utopia::stk::FunctionSpace;
    using FunctionSpaceTo = utopia::moonolith::FunctionSpace;

    auto params = param_list(param("mesh", param_list(param("type", "cube"), param("nx", nx))));
    FunctionSpaceFrom space_from;
    space_from.read(params);

    FunctionSpaceTo space_to;
    extract_trace_space(space_from, space_to);
    // utopia_test_assert(space_to.n_dofs() == (ny + 1) * (nz + 1) * 2);
}

void stk_moonolith_fe_transfer() {
    using FunctionSpace_t = utopia::stk::FunctionSpace;

    FunctionSpace_t space_from, space_to;

    auto params = param_list(param("mesh", param_list(param("type", "cube"))));

    space_from.read(params);
    space_to.read(params);

    utopia::stk::FETransfer transfer;
    utopia_test_assert(transfer.init(make_ref(space_from), make_ref(space_to)));
}

void interop_moonolith_stk() {
    UTOPIA_RUN_TEST(stk_moonolith_convert_mesh);
    UTOPIA_RUN_TEST(stk_moonolith_convert_space);
    UTOPIA_RUN_TEST(stk_moonolith_fe_transfer);
    UTOPIA_RUN_TEST(stk_moonolith_extract_surface);
    UTOPIA_RUN_TEST(stk_moonolith_extract_trace_space);
}

UTOPIA_REGISTER_TEST_FUNCTION(interop_moonolith_stk);
