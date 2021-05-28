#include "utopia_Testing.hpp"

#include "utopia_mars.hpp"
#include "utopia_ui.hpp"

using namespace utopia;

void mars_poisson() {
    using Mesh_t = utopia::mars::Mesh;
    using FunctionSpace_t = utopia::mars::FunctionSpace;
    using Vector = Traits<FunctionSpace_t>::Vector;
    using Matrix = Traits<FunctionSpace_t>::Matrix;

    Mesh_t mesh;
    mesh.unit_cube(2, 2, 2);

    FunctionSpace_t space;
    space.init(make_ref(mesh));

    space.add_dirichlet_boundary_condition("top", -0.1);
    space.add_dirichlet_boundary_condition("bottom", 0.1);

    Vector x, rhs;
    space.create_vector(x);
    space.create_vector(rhs);

    x.set(0.0);

    Matrix mat;
    space.create_matrix(mat);

    auto params = param_list(param("type", "LaplaceOperator"));

    OmniAssembler<FunctionSpace_t> assembler(make_ref(space));
    assembler.read(params);

    utopia_test_assert(assembler.assemble(x, mat, rhs));

    ConjugateGradient<Matrix, Vector> cg;

    utopia_test_assert(cg.solve(mat, rhs, x));
    utopia_test_assert(space.write("result.e", x));
}

void mars_assembler() { UTOPIA_RUN_TEST(mars_poisson); }

UTOPIA_REGISTER_TEST_FUNCTION(mars_assembler);
