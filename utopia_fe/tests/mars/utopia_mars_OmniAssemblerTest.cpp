#include "utopia_Testing.hpp"

#include "utopia_mars.hpp"
#include "utopia_ui.hpp"

#include "utopia_kokkos_LaplaceOperator_new.hpp"
#include "utopia_kokkos_UniformFE.hpp"
#include "utopia_mars_Discretization.hpp"
#include "utopia_mars_Material.hpp"

using namespace utopia;

using Mesh_t = utopia::mars::Mesh;
using FunctionSpace_t = utopia::mars::FunctionSpace;
using Vector_t = Traits<FunctionSpace_t>::Vector;
using Matrix_t = Traits<FunctionSpace_t>::Matrix;
using Scalar_t = Traits<FunctionSpace_t>::Scalar;

void mars_solve_aux(Input &in) {
    Chrono c;
    Chrono cc;
    c.start();
    cc.start();

    FunctionSpace_t space;
    space.read(in);

    cc.stop();
    utopia::out() << "Init: " << cc << "\n";

    int n_var = space.n_var();

    for (int v = 0; v < n_var; ++v) {
        space.add_dirichlet_boundary_condition("top", 0.1, v);
        space.add_dirichlet_boundary_condition("bottom", -0.1, v);
    }

    cc.start();
    Vector_t x, rhs;
    space.create_vector(x);
    space.create_vector(rhs);

    x.set(0.0);
    rhs.set(0.0);

    Matrix_t mat;
    space.create_matrix(mat);

    cc.stop();
    utopia::out() << "Create Matrix: " << cc << "\n";

    cc.start();
    OmniAssembler<FunctionSpace_t> assembler(make_ref(space));
    assembler.read(in);
    cc.stop();
    utopia::out() << "Prepare Assembly: " << cc << "\n";

    cc.start();
    utopia_test_assert(assembler.assemble(x, mat, rhs));

    cc.stop();
    utopia::out() << "Assembly: " << cc << "\n";

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
    cg.max_it(200);
    cg.verbose(true);

    utopia_test_assert(cg.solve(mat, rhs, x));

    c.stop();

    utopia::out() << "Solve: " << c << "\n";

    c.start();

    Scalar_t min_x = min(x);
    Scalar_t max_x = max(x);

    utopia::out() << "min_x: " << min_x << ", max_x: " << max_x << "\n";

    // static int once = 0;
    // if (once++ == 1) {
    //     utopia_test_assert(space.write("result.bp", x));
    // }

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

void mars_poisson_aux(int nx, int ny, int nz) {
    auto params =
        param_list(param("n_var", 1),
                   param("mesh", param_list(param("type", "cube"), param("nx", nx), param("ny", ny), param("nz", nz))),
                   param("material", param_list(param("type", "LaplaceOperator"))));

    mars_solve_aux(params);
}

void mars_linear_elasticity() { mars_linear_elasticity_aux(20, 20, 20); }

void mars_poisson_2D() { mars_poisson_aux(10, 11, 0); }

void mars_poisson_3D() { mars_poisson_aux(30, 24, 20); }

void mars_assembler() {
    UTOPIA_RUN_TEST(mars_poisson_2D);
    UTOPIA_RUN_TEST(mars_poisson_3D);

    // FIXME output.bp not working for vector values
    UTOPIA_RUN_TEST(mars_linear_elasticity);
}

UTOPIA_REGISTER_TEST_FUNCTION(mars_assembler);

using FS_t = utopia::mars::FunctionSpace;
using FE_t = utopia::kokkos::UniformFE<double>;
using Matrix_t = Traits<FS_t>::Matrix;
using Vector_t = Traits<FS_t>::Vector;
using Scalar_t = Traits<FS_t>::Scalar;
// using Assembler_t = utopia::mars::FEAssemblerNew<FE_t>;
using Discretization_t = utopia::Discretization<FS_t, FE_t>;
using Solver_t = utopia::ConjugateGradient<Matrix_t, Vector_t, HOMEMADE>;

// std::shared_ptr<utopia::mars::FEAssemblerNew<FE_t>> make_assembler(
//     const std::shared_ptr<utopia::mars::FunctionSpace> &space,
//     const int order) {
//     return std::make_shared<utopia::mars::FEAssemblerNew<FE_t>>();
// }

void mars_new_assembler_test() {
    auto params =
        param_list(param("n_var", 1),
                   param("mesh", param_list(param("type", "cube"), param("nx", 4), param("ny", 4), param("nz", 4))),
                   param("material", param_list(param("type", "LaplaceOperator"))));

    FS_t space;
    space.read(params);

    // auto assembler = make_assembler(make_ref(space), 2);

    utopia::kokkos::LaplaceOperatorNew<FS_t, FE_t> lapl;
    lapl.initialize(make_ref(space));

    Matrix_t mat;
    space.create_matrix(mat);

    Vector_t x, g;
    space.create_vector(x);
    space.create_vector(g);

    x.set(0.0);

    utopia_test_assert(lapl.hessian(x, mat));
    utopia_test_assert(lapl.gradient(x, g));

    Scalar_t ng = norm2(g);
    Scalar_t nx = norm2(x);
    Scalar_t nm = norm2(mat);

    utopia_test_assert(nm > 0.0);

    g *= -1;
    space.apply_constraints(mat, g);
    space.apply_constraints(x);

    Solver_t solver;
    solver.apply_gradient_descent_step(true);
    utopia_test_assert(solver.solve(mat, g, x));
}

UTOPIA_REGISTER_TEST_FUNCTION(mars_new_assembler_test);
