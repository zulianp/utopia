#include "utopia_fe_base.hpp"

#ifdef UTOPIA_WITH_INTREPID2

#include "utopia_Testing.hpp"
#include "utopia_ui.hpp"

#include "utopia_FEInteroperability.hpp"

// Stk includes
#include "utopia_libmesh_FunctionSpace_new.hpp"
#include "utopia_libmesh_Mesh.hpp"

// Intrepid2 includes
// #include "utopia_intrepid2_LaplaceOperator.hpp"
// #include "utopia_intrepid2_LinearElasticity.hpp"
// #include "utopia_intrepid2_NeoHookean.hpp"
// #include "utopia_intrepid2_ShellTools.hpp"
// #include "utopia_intrepid2_VectorLaplaceOperator.hpp"

#include "utopia_libmesh_kokkos.hpp"

// FIXME: This is the last include because the operator files are not yet in the correct place
#include "utopia_SpaceAndFETest.hpp"

using namespace utopia;

using LibMeshScalar_t = Traits<utopia::libmesh::FunctionSpace>::Scalar;
void interop_libmesh_kokkos() {
    SpaceAndFETest<utopia::libmesh::FunctionSpace, utopia::kokkos::FE<LibMeshScalar_t>>().run();
}

void interop_libmesh_kokkos_test_create_fe_boundary() {
    utopia::libmesh::FunctionSpace space;
    // auto params = param_list(
    //     param("mesh", param_list(param("nx", 3), param("ny", 3), param("nz", 3), param("elem_type", "HEX8"))));

    auto params = param_list(param(
        "mesh", param_list(param("type", "square"), param("nx", 3), param("ny", 3), param("elem_type", "QUAD4"))));

    space.read(params);

    space.describe(utopia::out().stream());

    utopia::kokkos::FE<LibMeshScalar_t> fe;
    create_fe(space, fe, 2);
    // create_fe_on_boundary(space, fe, 2);
    // create_fe_on_boundary(space, fe, "top", 2);

    // std::cout << "fe.n_cells(): " << fe.n_cells() << std::endl;

    // space.mesh().write("mesh.e");

    LibMeshViewDevice_t evec("evec", fe.n_cells(), fe.n_shape_functions());
    Kokkos::deep_copy(evec, 1.0);

    PetscVector vec;
    space.create_vector(vec);
    local_to_global(space, evec, utopia::ADD_MODE, vec);

    double norm_vec = norm1(vec);
    disp(norm_vec);
}

UTOPIA_REGISTER_TEST_FUNCTION(interop_libmesh_kokkos);
UTOPIA_REGISTER_TEST_FUNCTION(interop_libmesh_kokkos_test_create_fe_boundary);

#endif  // UTOPIA_WITH_INTREPID2
