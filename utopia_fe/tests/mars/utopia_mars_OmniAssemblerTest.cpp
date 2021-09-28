#include "utopia_Testing.hpp"

#include "utopia_mars.hpp"
#include "utopia_ui.hpp"

using namespace utopia;

using Mesh_t = utopia::mars::Mesh;
using FunctionSpace_t = utopia::mars::FunctionSpace;
using Vector_t = Traits<FunctionSpace_t>::Vector;
using Matrix_t = Traits<FunctionSpace_t>::Matrix;
using Scalar_t = Traits<FunctionSpace_t>::Scalar;

void mars_solve_aux(Input &in) {
    Chrono c;
    c.start();

    FunctionSpace_t space;
    space.read(in);
    int n_var = space.n_var();

    for (int v = 0; v < n_var; ++v) {
        space.add_dirichlet_boundary_condition("top", 0.1, v);
        space.add_dirichlet_boundary_condition("bottom", -0.1, v);
    }

    Vector_t x, rhs;
    space.create_vector(x);
    space.create_vector(rhs);

    x.set(0.0);
    rhs.set(0.0);

    Matrix_t mat;
    space.create_matrix(mat);

    OmniAssembler<FunctionSpace_t> assembler(make_ref(space));
    assembler.read(in);

    utopia_test_assert(assembler.assemble(x, mat, rhs));

    Vector_t ones(layout(x), 1.0);
    Vector_t sum_rows = mat * ones;
    Scalar_t sum_mat = sum(abs(sum_rows));

    // write("load_mat.mm", mat);
    // write("load_rhs.mm", rhs);

    utopia_test_assert(sum_mat < 1e-8);

    space.apply_constraints(mat, rhs);
    // space.apply_constraints(x);

    c.stop();

    utopia::out() << "Discretization: " << c << "\n";

    c.start();

    ConjugateGradient<Matrix_t, Vector_t> cg;
    cg.set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix_t, Vector_t>>());
    cg.apply_gradient_descent_step(true);
    cg.verbose(true);

    utopia_test_assert(cg.solve(mat, rhs, x));

    c.stop();

    utopia::out() << "Solve: " << c << "\n";

    c.start();

    Scalar_t min_x = min(x);
    Scalar_t max_x = max(x);

    utopia::out() << "min_x: " << min_x << ", max_x: " << max_x << "\n";

    // x.set(1.);

    utopia_test_assert(
        space.write("result.vtu",
                    // "result",
                    // + std::to_string(nx) + "x" + std::to_string(ny) + "x" + std::to_string(nz) + ".vtu",
                    x));

    c.stop();

    utopia::out() << "Output: " << c << "\n";
}

void mars_linear_elasticity_aux(int nx, int ny, int nz) {
    auto params =
        param_list(param("n_var", 2 + (nz != 0)),
                   param("mesh", param_list(param("type", "cube"), param("nx", nx), param("ny", ny), param("nz", nz))),
                   param("material", param_list(param("type", "LinearElasticity"))));

    mars_solve_aux(params);
}

void mars_linear_elasticity() { mars_linear_elasticity_aux(40, 40, 40); }

void mars_poisson_aux(int nx, int ny, int nz) {
    using Mesh_t = utopia::mars::Mesh;
    using FunctionSpace_t = utopia::mars::FunctionSpace;
    using Vector_t = Traits<FunctionSpace_t>::Vector;
    using Matrix_t = Traits<FunctionSpace_t>::Matrix;
    using Scalar_t = Traits<FunctionSpace_t>::Scalar;

    Chrono c;
    c.start();

    Mesh_t mesh;
    mesh.unit_cube(nx, ny, nz);

    FunctionSpace_t space;
    InputParameters space_params;
    space_params.set("n_var", 1);
    space.read(space_params);
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

    auto params = param_list(param("material", param_list(param("type", "LaplaceOperator"))));

    OmniAssembler<FunctionSpace_t> assembler(make_ref(space));
    assembler.read(params);

    utopia_test_assert(assembler.assemble(x, mat, rhs));

    Vector_t ones(layout(x), 1.0);
    Vector_t sum_rows = mat * ones;
    Scalar_t sum_mat = sum(abs(sum_rows));

    write("load_mat.mm", mat);
    write("load_rhs.mm", rhs);

    utopia_test_assert(sum_mat < 1e-8);

    space.apply_constraints(mat, rhs);
    space.apply_constraints(x);

    c.stop();

    utopia::out() << "Discretization: " << c << "\n";

    c.start();

    ConjugateGradient<Matrix_t, Vector_t> cg;
    cg.set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix_t, Vector_t>>());
    cg.apply_gradient_descent_step(true);
    cg.verbose(true);

    utopia_test_assert(cg.solve(mat, rhs, x));

    c.stop();

    utopia::out() << "Solve: " << c << "\n";

    c.start();

    utopia_test_assert(
        space.write("result" + std::to_string(nx) + "x" + std::to_string(ny) + "x" + std::to_string(nz) + ".vtu", x));

    c.stop();

    utopia::out() << "Output: " << c << "\n";
}

void mars_poisson_2D() { mars_poisson_aux(10, 11, 0); }

void mars_poisson_3D() { mars_poisson_aux(8, 9, 10); }

void mars_assembler() {
    UTOPIA_RUN_TEST(mars_poisson_2D);
    UTOPIA_RUN_TEST(mars_poisson_3D);

    UTOPIA_RUN_TEST(mars_linear_elasticity);
}

UTOPIA_REGISTER_TEST_FUNCTION(mars_assembler);
