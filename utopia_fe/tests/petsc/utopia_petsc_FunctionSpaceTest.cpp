#include "utopia_Testing.hpp"
#include "utopia_UnitTest.hpp"

#include "utopia_petsc_FunctionSpace.hpp"
#include "utopia_ui.hpp"

using namespace utopia;

void petsc_create_matrix() {
    using Mesh_t = utopia::petsc::StructuredGrid;
    using FunctionSpace_t = utopia::petsc::FunctionSpace;
    using Communicator = Traits<FunctionSpace_t>::Communicator;
    using Matrix = Traits<FunctionSpace_t>::Matrix;

    Communicator comm;

    Mesh_t mesh(comm);
    mesh.unit_cube(2, 2, 2);

    FunctionSpace_t space(make_ref(mesh));

    Matrix mat;
    space.create_matrix(mat);

    SizeType n = mat.rows();

    utopia_test_assert(!mat.empty());
    UTOPIA_TEST_EQ(n, mesh.n_nodes());
}

void petsc_create_vector() {
    using Mesh_t = utopia::petsc::StructuredGrid;
    using FunctionSpace_t = utopia::petsc::FunctionSpace;
    using Communicator = Traits<FunctionSpace_t>::Communicator;
    using Vector = Traits<FunctionSpace_t>::Vector;
    using SizeType = Traits<FunctionSpace_t>::SizeType;

    Communicator comm;

    Mesh_t mesh(comm);
    mesh.unit_cube(2, 2, 2);
    FunctionSpace_t space(make_ref(mesh));

    Vector vec;
    space.create_vector(vec);

    SizeType n = vec.size();

    utopia_test_assert(!vec.empty());
    UTOPIA_TEST_EQ(n, mesh.n_nodes());
}

void petsc_function_space() {
    UTOPIA_RUN_TEST(petsc_create_matrix);
    UTOPIA_RUN_TEST(petsc_create_vector);
}

UTOPIA_REGISTER_TEST_FUNCTION(petsc_function_space);
