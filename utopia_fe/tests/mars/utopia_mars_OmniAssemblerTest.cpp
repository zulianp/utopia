#include "utopia_Testing.hpp"

#include "utopia_mars.hpp"
#include "utopia_ui.hpp"

using namespace utopia;

void mars_assemble_poisson() {
    using Mesh_t = utopia::mars::Mesh;
    using FunctionSpace_t = utopia::mars::FunctionSpace;
    using Vector = Traits<FunctionSpace_t>::Vector;
    using Matrix = Traits<FunctionSpace_t>::Matrix;

    Mesh_t mesh;
    mesh.unit_cube(2, 2, 2);

    FunctionSpace_t space;
    space.init(make_ref(mesh));

    Vector x, rhs;
    space.create_vector(x);
    space.create_vector(rhs);

    Matrix mat;
    space.create_matrix(mat);

    OmniAssembler<FunctionSpace_t> assembler(make_ref(space));
    utopia_test_assert(assembler.assemble(x, mat, rhs));
}

void mars_assembler() { UTOPIA_RUN_TEST(mars_assemble_poisson); }

UTOPIA_REGISTER_TEST_FUNCTION(mars_assembler);
