#include "utopia_Testing.hpp"

#include "utopia_mars.hpp"
#include "utopia_ui.hpp"

using namespace utopia;

void mars_create_matrix() {
    using Mesh_t = utopia::mars::Mesh;
    using FunctionSpace_t = utopia::mars::FunctionSpace;
    using Matrix = Traits<FunctionSpace_t>::Matrix;

    Mesh_t mesh;
    mesh.unit_cube(2, 2, 2);

    FunctionSpace_t space;
    space.init(make_ref(mesh));

    Matrix mat;
    space.create_matrix(mat);

    // space.describe();

    utopia_test_assert(!mat.empty());
}

void mars_create_vector() {
    using Mesh_t = utopia::mars::Mesh;
    using FunctionSpace_t = utopia::mars::FunctionSpace;
    using Vector = Traits<FunctionSpace_t>::Vector;

    Mesh_t mesh;
    mesh.unit_cube(2, 2, 2);

    FunctionSpace_t space;
    space.init(make_ref(mesh));

    Vector vec;
    space.create_vector(vec);

    utopia_test_assert(!vec.empty());
}

void mars_function_space() {
    UTOPIA_RUN_TEST(mars_create_matrix);
    UTOPIA_RUN_TEST(mars_create_vector);
}

UTOPIA_REGISTER_TEST_FUNCTION(mars_function_space);
