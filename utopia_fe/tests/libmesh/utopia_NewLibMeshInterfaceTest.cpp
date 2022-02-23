#include "utopia_Testing.hpp"

#include "utopia_MeshTest.hpp"
#include "utopia_RunParallelTest.hpp"

#include "utopia_libmesh_FunctionSpace_new.hpp"
#include "utopia_libmesh_Mesh.hpp"
#include "utopia_ui.hpp"

using namespace utopia;

void lm_create_mesh() {
    using Mesh_t = utopia::libmesh::Mesh;

    InputParameters params;
    params.set("type", "cube");
    params.set("elem_type", "TET4");

    Mesh_t mesh;
    mesh.read(params);
}

void lm_create_function_space() {
    using FunctionSpace_t = utopia::libmesh::FunctionSpace;

    auto params = param_list(                     //
        param("system_type", "linear_implicit"),  //
        param("mesh",                             // Nested level of settings
              param_list(param("type", "sphere"))));

    FunctionSpace_t space;
    space.read(params);

    std::cout << space.n_local_dofs() << " " << space.n_dofs() << std::endl;

    UVector v;
    USparseMatrix m;

    space.create_vector(v);
    space.create_matrix(m);
}

void lm_read_function_space_with_data() {
    using FunctionSpace_t = utopia::libmesh::FunctionSpace;

    FunctionSpace_t space;

    UVector displacement;
    space.read("squeezed_pump.e", {"disp_x", "disp_y", "disp_z"}, displacement);
    space.mesh().describe(std::cout);

    double norm_displacement = norm2(displacement);
    disp(norm_displacement);

    // std::cout << space.n_local_dofs() << " " << space.n_dofs() << std::endl;

    // space.mesh().write("lm_read_function_space_with_data.e");
    space.write("squeezed_pump_rw.e", displacement);
}

void lm() {
    UTOPIA_RUN_TEST(lm_create_mesh);
    UTOPIA_RUN_TEST(lm_create_function_space);
    // UTOPIA_RUN_TEST(lm_read_function_space_with_data);

    utopia::run_parallel_test<MeshTest<utopia::libmesh::Mesh>>();
}

UTOPIA_REGISTER_TEST_FUNCTION(lm);
