#include "utopia_Testing.hpp"

#include "utopia_mars.hpp"
#include "utopia_ui.hpp"

using namespace utopia;

void mars_poisson() {
    using Mesh_t = utopia::mars::Mesh;
    using FunctionSpace_t = utopia::mars::FunctionSpace;
    using Vector_t = Traits<FunctionSpace_t>::Vector;
    using Matrix_t = Traits<FunctionSpace_t>::Matrix;
    using Scalar_t = Traits<FunctionSpace_t>::Scalar;

    int n = 400;
    Mesh_t mesh;
    mesh.unit_cube(n, n, 0);

    FunctionSpace_t space;
    space.init(make_ref(mesh));

    space.add_dirichlet_boundary_condition("top", -0.1);
    space.add_dirichlet_boundary_condition("bottom", 0.1);

    Vector_t x, rhs;
    space.create_vector(x);
    space.create_vector(rhs);

    x.set(0.0);
    rhs.set(0.0);

    Matrix_t mat;
    space.create_matrix(mat);

    auto params = param_list(param("type", "LaplaceOperator"));

    OmniAssembler<FunctionSpace_t> assembler(make_ref(space));
    assembler.read(params);

    utopia_test_assert(assembler.assemble(x, mat, rhs));

    Vector_t ones(layout(x), 1.0);
    Vector_t sum_rows = mat * ones;
    Scalar_t sum_mat = sum(abs(sum_rows));

    utopia_test_assert(sum_mat < 1e-8);

    // write("load_mat.mm", mat);
    // write("load_rhs.mm", rhs);

    space.apply_constraints(mat, rhs);
    // space.apply_constraints(x);

    ConjugateGradient<Matrix_t, Vector_t> cg;
    cg.set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix_t, Vector_t>>());
    cg.apply_gradient_descent_step(true);
    cg.verbose(true);

    utopia_test_assert(cg.solve(mat, rhs, x));
    utopia_test_assert(space.write("result.vtu", x));
}

void mars_assembler() { UTOPIA_RUN_TEST(mars_poisson); }

UTOPIA_REGISTER_TEST_FUNCTION(mars_assembler);
