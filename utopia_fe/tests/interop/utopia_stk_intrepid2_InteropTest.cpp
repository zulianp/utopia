#include "utopia_fe_base.hpp"

#ifdef UTOPIA_WITH_INTREPID2

#include "utopia_Testing.hpp"
#include "utopia_ui.hpp"

#include "utopia_stk_FunctionSpace.hpp"
#include "utopia_stk_Mesh.hpp"

#include "utopia_FEInteroperability.hpp"

#include "utopia_intrepid2_LaplaceOperator.hpp"

using namespace utopia;

void poisson_problem() {
    using FunctionSpace_t = utopia::stk::FunctionSpace;
    using Scalar_t = Traits<FunctionSpace_t>::Scalar;
    using Matrix_t = Traits<FunctionSpace_t>::Matrix;
    using Vector_t = Traits<FunctionSpace_t>::Vector;
    using FE_t = utopia::intrepid2::FE<Scalar_t>;

    auto params =
        param_list(param("path",
                         // "../data/knf/pump/membrane.e"
                         // "../data//bugs/transfer/cube_10x10x10.e"
                         "/Users/zulianp/Desktop/code/fluyafsi/build_opt/02_Coarser_Thick/res_pitzDaily_coupled.e"
                         // "../data/knf/rectangle_4_tris.e"
                         ));

    FunctionSpace_t space;
    space.read(params);
    // space.add_dirichlet_boundary_condition("surface_1", 1.0);
    // space.add_dirichlet_boundary_condition("surface_3", -1.0);

    space.add_dirichlet_boundary_condition("inlet", 1.0);
    space.add_dirichlet_boundary_condition("outlet", -1.0);

    // space.describe(std::cout);

    // Chrono c;
    // c.start();

    auto fe_ptr = std::make_shared<FE_t>();
    create_fe(space, *fe_ptr, 2);

    LaplaceOperator<Scalar_t> lapl{1.0};
    utopia::intrepid2::Assemble<LaplaceOperator<Scalar_t>> assembler(lapl, fe_ptr);
    assembler.init();

    // local to global
    Matrix_t mat;
    local_to_global(space, assembler.element_matrices(), mat);

    Vector_t rhs, x;
    space.create_vector(rhs);

    x.zeros(layout(rhs));

    space.apply_constraints(mat, rhs);

    Factorization<Matrix_t, Vector_t> solver;

    utopia_test_assert(solver.solve(mat, rhs, x));

    space.write("poisson_problem.e", x);

    write("poisson_problem.m", mat);
    write("poisson_problem_rhs.m", rhs);
    write("poisson_problem_x.m", x);

    // c.stop();
    // std::cout << c << std::endl;

    // std::ofstream os("prova.txt");
    // assembler.describe(os);
    // os.close();
}

void interop_stk_intrepid2() { UTOPIA_RUN_TEST(poisson_problem); }

UTOPIA_REGISTER_TEST_FUNCTION(interop_stk_intrepid2);

#endif  // UTOPIA_WITH_INTREPID2