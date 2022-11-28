#include "utopia_fe_base.hpp"

#ifdef UTOPIA_WITH_INTREPID2

#include "utopia_Testing.hpp"
#include "utopia_ui.hpp"

#include "utopia_FEInteroperability.hpp"

// Stk includes
#include "utopia_stk_FunctionSpace.hpp"
#include "utopia_stk_Mesh.hpp"

// utopia/kokkos includes
#include "utopia_kokkos_L2Norm.hpp"
#include "utopia_kokkos_LaplaceOperator.hpp"
#include "utopia_kokkos_LinearElasticity.hpp"
#include "utopia_kokkos_NeoHookean.hpp"
#include "utopia_kokkos_VectorLaplaceOperator.hpp"

#include "utopia_hyperelasticity_NeoHookeanOgden_3.hpp"
#include "utopia_kokkos_AutoHyperElasticityNew.hpp"

// utopia/intrepid2 includes
#include "utopia_intrepid2_ShellTools.hpp"

// Interop includes
#include "utopia_stk_intrepid2.hpp"

// FIXME: This is the last include because the operator files are not yet in the correct place
#include "utopia_SpaceAndFETest.hpp"

#include "utopia_ScalarProductTest.hpp"

#include "utopia_kokkos_LaplaceOperator_new.hpp"

#include "utopia_stk_intrepid2_Discretization.hpp"
#include "utopia_stk_intrepid2_Material.hpp"

using namespace utopia;

using StkScalar = Traits<utopia::stk::FunctionSpace>::Scalar;
void interop_stk_intrepid2() { SpaceAndFETest<utopia::stk::FunctionSpace, utopia::intrepid2::FE<StkScalar>>().run(); }

UTOPIA_REGISTER_TEST_FUNCTION(interop_stk_intrepid2);

void scalar_product_stk_intrepid2() { ScalarProductTest<utopia::stk::FunctionSpace>().run(); }

UTOPIA_REGISTER_TEST_FUNCTION(scalar_product_stk_intrepid2);

void new_assembler_test() {
    using FS_t = utopia::stk::FunctionSpace;
    using FE_t = utopia::intrepid2::FE<double>;
    using Matrix_t = Traits<FS_t>::Matrix;
    using Vector_t = Traits<FS_t>::Vector;
    using Scalar_t = Traits<FS_t>::Scalar;
    using Assembler_t = utopia::kokkos::FEAssembler<FS_t, FE_t>;
    using Discretization_t = utopia::Discretization<FS_t, FE_t>;
    using Solver_t = utopia::ConjugateGradient<Matrix_t, Vector_t, HOMEMADE>;

    UnitCubeSpaceAndFETest<FS_t, FE_t> test;
    auto params = test.cube_space_param(1);

    FS_t space;
    space.read(params);
    test.add_cube_bc(space, 1);

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
    solver.verbose(true);
    utopia_test_assert(solver.solve(mat, g, x));
}

UTOPIA_REGISTER_TEST_FUNCTION(new_assembler_test);

void new_auto_assembler_test() {
    using FS_t = utopia::stk::FunctionSpace;
    using FE_t = utopia::intrepid2::FE<double>;
    using Matrix_t = Traits<FS_t>::Matrix;
    using Vector_t = Traits<FS_t>::Vector;
    using Scalar_t = Traits<FS_t>::Scalar;
    using Assembler_t = utopia::kokkos::FEAssembler<FS_t, FE_t>;
    using Discretization_t = utopia::Discretization<FS_t, FE_t>;
    using Solver_t = utopia::ConjugateGradient<Matrix_t, Vector_t, HOMEMADE>;

    int n = 80;
    auto params =
        param_list(param("n_var", 3),
                   param("mesh", param_list(param("type", "cube"), param("nx", n), param("ny", n), param("nz", n))),
                   param("material", param_list(param("type", "LaplaceOperator"))));

    FS_t space;
    space.read(params);

    auto l = utopia::stk::SideSet::Cube::convert("left");
    auto r = utopia::stk::SideSet::Cube::convert("right");

    space.add_dirichlet_boundary_condition(l, -0.05, 0);
    space.add_dirichlet_boundary_condition(r, 0.05, 0);

    space.add_dirichlet_boundary_condition(l, 0, 1);
    space.add_dirichlet_boundary_condition(r, 0.05, 1);

    space.add_dirichlet_boundary_condition(l, 0, 2);
    space.add_dirichlet_boundary_condition(r, 0, 2);

    utopia::kokkos::AutoHyperElasticityNew<FS_t, FE_t, utopia::kernels::NeoHookeanOgden<Scalar_t, 3>> neohook;
    neohook.initialize(make_ref(space));

    Matrix_t mat;
    space.create_matrix(mat);

    Vector_t x, g;
    space.create_vector(x);
    space.create_vector(g);

    x.set(0.0);

    utopia_test_assert(neohook.hessian(x, mat));
    utopia_test_assert(neohook.gradient(x, g));

    Scalar_t ng = norm2(g);
    Scalar_t nx = norm2(x);
    Scalar_t nm = norm2(mat);

    Scalar_t sm = sum(mat);
    Scalar_t sg = sum(g);

    utopia_test_assert(sg < 1e-8);
    utopia_test_assert(sm < 1e-8);
    utopia_test_assert(nm > 0.0);

    g *= -1;
    space.apply_constraints(mat, g);
    space.apply_constraints(x);

    Solver_t solver;
    solver.apply_gradient_descent_step(true);
    solver.verbose(true);
    utopia_test_assert(solver.solve(mat, g, x));

    space.write("neo.e", x);
}

UTOPIA_REGISTER_TEST_FUNCTION(new_auto_assembler_test);

#endif  // UTOPIA_WITH_INTREPID2
