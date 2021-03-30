#include "utopia_fe_base.hpp"

#ifdef UTOPIA_WITH_INTREPID2

#include "utopia_Testing.hpp"
#include "utopia_ui.hpp"

#include "utopia_stk_FunctionSpace.hpp"
#include "utopia_stk_Mesh.hpp"

#include "utopia_FEInteroperability.hpp"

#include "utopia_intrepid2_LaplaceOperator.hpp"
#include "utopia_intrepid2_LinearElasticity.hpp"
#include "utopia_intrepid2_VectorLaplaceOperator.hpp"

using namespace utopia;

class FEProblemTest {
public:
    using FunctionSpace_t = utopia::stk::FunctionSpace;
    using Scalar_t = Traits<FunctionSpace_t>::Scalar;
    using Matrix_t = Traits<FunctionSpace_t>::Matrix;
    using Vector_t = Traits<FunctionSpace_t>::Vector;
    using FE_t = utopia::intrepid2::FE<Scalar_t>;

    bool save_output{false};
    bool export_tensors{false};
    bool verbose{false};
    Scalar_t rtol{1e-6};

    template <class Op>
    void assemble_and_solve(const std::string &name, FunctionSpace_t &space, Op &op) {
        auto fe_ptr = std::make_shared<FE_t>();
        create_fe(space, *fe_ptr, 2);

        utopia::intrepid2::Assemble<Op> assembler(op, fe_ptr);
        assembler.init();

        Matrix_t mat;
        local_to_global(space, assembler.element_matrices(), mat);

        Vector_t row_sum = sum(mat, 1);
        Scalar_t sum_row_sum = sum(abs(row_sum));
        utopia_test_assert(approxeq(sum_row_sum, 0.0, 1e-9));

        Vector_t rhs, x;
        space.create_vector(rhs);

        x.zeros(layout(rhs));

        space.apply_constraints(mat, rhs);

        ////////////////////////////////////////////////////////////

        if (export_tensors) {
            write(name + ".m", mat);
            write(name + "_rhs.m", rhs);

            std::ofstream os(name + ".txt");
            assembler.describe(os);
            os.close();
        }

        ////////////////////////////////////////////////////////////

        KSPSolver<Matrix_t, Vector_t> solver;
        solver.pc_type("hypre");
        // solver.pc_type(PCILU);

        // if (mat.is_block()) {
        //     solver.factor_set_pivot_in_blocks(true);
        // }

        solver.rtol(rtol);
        solver.verbose(verbose);

        utopia_test_assert(solver.solve(mat, rhs, x));

        ////////////////////////////////////////////////////////////

        if (save_output) {
            space.write(name + ".e", x);
            // space.write(name + ".e", rhs);
        }

        if (export_tensors) {
            write(name + "_x.m", x);
        }
    }

    void poisson_problem() {
        auto params = param_list(param(
            "mesh",
            param_list(param(
                "path", "/Users/zulianp/Desktop/code/fluyafsi/build_opt/02_Coarser_Thick/res_pitzDaily_coupled.e"))));

        FunctionSpace_t space;
        space.read(params);
        space.add_dirichlet_boundary_condition("inlet", 1.0);
        space.add_dirichlet_boundary_condition("outlet", -1.0);

        LaplaceOperator<Scalar_t> lapl{1.0};

        assemble_and_solve("poisson", space, lapl);
    }

    void vector_poisson_problem() {
        auto params = param_list(
            param("n_var", 3),
            param("mesh",
                  param_list(param(
                      "path",
                      "/Users/zulianp/Desktop/code/fluyafsi/build_opt/02_Coarser_Thick/res_pitzDaily_coupled.e"))));

        FunctionSpace_t space;
        space.read(params);
        space.add_dirichlet_boundary_condition("inlet", 1.0, 0);
        space.add_dirichlet_boundary_condition("inlet", 2.0, 1);
        space.add_dirichlet_boundary_condition("inlet", 3.0, 2);

        space.add_dirichlet_boundary_condition("outlet", -1.0, 0);
        space.add_dirichlet_boundary_condition("outlet", -2.0, 1);
        space.add_dirichlet_boundary_condition("outlet", -3.0, 2);

        VectorLaplaceOperator<3, Scalar_t> lapl{1.0};

        assemble_and_solve("vector_poisson", space, lapl);
    }

    void elasticity_problem() {
        static const int Dim = 3;
        auto params = param_list(
            param("n_var", Dim),
            param("mesh",
                  param_list(param(
                      "path",
                      "/Users/zulianp/Desktop/code/fluyafsi/build_opt/02_Coarser_Thick/res_pitzDaily_coupled.e"))));

        FunctionSpace_t space;
        space.read(params);
        space.add_dirichlet_boundary_condition("inlet", -0.005, 0);
        space.add_dirichlet_boundary_condition("inlet", 0.0, 1);
        space.add_dirichlet_boundary_condition("inlet", 0.001, 2);

        space.add_dirichlet_boundary_condition("outlet", 0.005, 0);
        space.add_dirichlet_boundary_condition("outlet", 0.0, 1);
        space.add_dirichlet_boundary_condition("outlet", -0.001, 2);

        LinearElasticity<Dim, Scalar_t> linear_elasticity{1.0, 1.0};
        assemble_and_solve("elasticity", space, linear_elasticity);
    }

    void poisson_problem_parallel() {
        // auto params = param_list(param("mesh", param_list(param("path", "../data/knf/cube_vs_cube/body.e"))));
        auto params = param_list(param("mesh", param_list(param("path", "../data/knf/rectangle_4_tris.e"))));

        FunctionSpace_t space;
        space.read(params);
        space.add_dirichlet_boundary_condition("body_top", 1.0);
        space.add_dirichlet_boundary_condition("body_bottom", -1.0);

        LaplaceOperator<Scalar_t> lapl{1.0};

        std::stringstream ss;
        space.describe(ss);

        // if (space.comm().size() == 2) {
        //     space.comm().synched_print(ss.str());
        // }

        assemble_and_solve("poisson_problem_parallel", space, lapl);
    }

    void elasticity_problem_parallel() {
        static const int Dim = 3;
        auto params = param_list(param("n_var", Dim),
                                 param("mesh", param_list(param("path", "../data/knf/cube_vs_cube/body.e"))));

        FunctionSpace_t space;
        space.read(params);
        space.add_dirichlet_boundary_condition("body_top", 0.0, 0);
        space.add_dirichlet_boundary_condition("body_top", -0.1, 1);
        space.add_dirichlet_boundary_condition("body_top", 0.0, 2);

        space.add_dirichlet_boundary_condition("body_bottom", 0.0, 0);
        space.add_dirichlet_boundary_condition("body_bottom", 0.1, 1);
        space.add_dirichlet_boundary_condition("body_bottom", 0.0, 2);

        LinearElasticity<Dim, Scalar_t> linear_elasticity{1.0, 1.0};
        assemble_and_solve("elasticity_problem_parallel", space, linear_elasticity);
    }

    void run() {
        if (mpi_world_size() == 1) {
            UTOPIA_RUN_TEST(poisson_problem);
            UTOPIA_RUN_TEST(vector_poisson_problem);
            UTOPIA_RUN_TEST(elasticity_problem);
        }

        if (mpi_world_size() <= 2) {
            save_output = true;
            UTOPIA_RUN_TEST(poisson_problem_parallel);
            // UTOPIA_RUN_TEST(elasticity_problem_parallel);
        }
    }
};

void interop_stk_intrepid2() { FEProblemTest().run(); }

UTOPIA_REGISTER_TEST_FUNCTION(interop_stk_intrepid2);

#endif  // UTOPIA_WITH_INTREPID2
