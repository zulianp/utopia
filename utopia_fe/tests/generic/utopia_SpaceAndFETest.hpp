#ifndef UTOPIA_SPACE_AND_FE_TEST_HPP
#define UTOPIA_SPACE_AND_FE_TEST_HPP

#include "utopia_Testing.hpp"
#include "utopia_fe_base.hpp"
#include "utopia_ui.hpp"

#include "utopia_fe_Core.hpp"

namespace utopia {

    template <class FunctionSpace, class FE>
    class SpaceAndFETest {
    public:
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using SideSet_t = typename Traits<FunctionSpace>::SideSet;
        using Field_t = utopia::Field<FunctionSpace>;

        bool save_output{false};
        bool export_tensors{false};
        bool verbose{false};
        bool use_cube{true};
        Scalar_t rtol{1e-6};

        static std::string get_shell_mesh_path() { return "../data/fe_problem_solve/shell.e"; }

        static std::string get_cube_path() {
            std::string dir = "../data/knf/cube_vs_cube";
            if (mpi_world_size() > 1) {
                dir += "/" + std::to_string(mpi_world_size());
            }

            return dir + "/body.e";
        }

        InputParameters cube_space_param(const int n_var) const {
            const int nx = 10;
            const int ny = 10;
            const int nz = 10;

            return param_list(
                param("n_var", n_var),
                param("mesh", param_list(param("type", "cube"), param("nx", nx), param("ny", ny), param("nz", nz))));
        }

        InputParameters input_params(const std::string &path, const int n_var = 1) const {
            if (use_cube) {
                return cube_space_param(n_var);
            } else {
                return param_list(param("n_var", 3),
                                  param("mesh", param_list(param("type", "file"), param("path", path))));
            }
        }

        inline void add_cube_bc(FunctionSpace &space, const int n_var = 1) const {
            for (int i = 0; i < n_var; ++i) {
                space.add_dirichlet_boundary_condition(SideSet_t::top(), (i == 0), i);
                space.add_dirichlet_boundary_condition(SideSet_t::bottom(), -(i == 0), i);
            }
        }

        static std::string get_more_complex_mesh_path() {
            std::string dir = "/Users/zulianp/Desktop/code/utopia/utopia_fe/data/knf/pitzdaily/";
            if (mpi_world_size() > 1) {
                dir += "/" + std::to_string(mpi_world_size());
            }
            return dir + "/pitz_daily.e";
        }

        static std::string get_2D_mesh_path() {
            std::string dir = "../data/knf/";
            if (mpi_world_size() > 1) {
                dir += "/" + std::to_string(mpi_world_size());
            }

            return dir + "/rectangle_4_tris.e";
        }

        template <class Op>
        void assemble_and_solve(const std::string &name, FunctionSpace &space, Op &op) {
            using Assembler_t = typename AssembleTraits<FE, Op>::Type;
            using FEField_t = typename Traits<FE>::Field;

            auto fe_ptr = std::make_shared<FE>();
            create_fe(space, *fe_ptr, 2);

            Vector_t rhs;
            space.create_vector(rhs);

            Field_t field_x;
            space.create_field(field_x);

            auto &x = field_x.data();
            x.set(0.0);

            FEField_t fe_field(fe_ptr);
            convert_field(field_x, fe_field);

            Assembler_t assembler(fe_ptr, op);
            assembler.update(make_ref(fe_field));
            assembler.assemble();

            Matrix_t mat;

            assert(assembler.is_matrix());

            local_to_global(space, assembler.matrix_data(), OVERWRITE_MODE, mat);

            if (assembler.is_vector()) {
                local_to_global(space, assembler.vector_data(), OVERWRITE_MODE, rhs);
            }

            Vector_t row_sum = sum(mat, 1);
            Scalar_t sum_row_sum = sum(abs(row_sum));

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

            utopia_test_assert(approxeq(sum_row_sum, 0.0, 1e-9));

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

        void create_fe_test() {
            auto params = cube_space_param(1);
            FunctionSpace space;
            space.read(params);

            auto fe_ptr = std::make_shared<FE>();
            create_fe(space, *fe_ptr, 0);
        }

        void poisson_problem() {
            auto params = input_params(get_more_complex_mesh_path());
            FunctionSpace space;
            space.read(params);

            if (use_cube) {
                add_cube_bc(space, 1);
            } else {
                space.add_dirichlet_boundary_condition("inlet", 1.0);
                space.add_dirichlet_boundary_condition("outlet", -1.0);
            }

            LaplaceOperator<Scalar_t> lapl{1.0};

            assemble_and_solve("poisson", space, lapl);
        }

        void vector_poisson_problem() {
            auto params = input_params(get_more_complex_mesh_path(), 3);

            FunctionSpace space;
            space.read(params);

            if (use_cube) {
                add_cube_bc(space, 3);
            } else {
                space.add_dirichlet_boundary_condition("inlet", 1.0, 0);
                space.add_dirichlet_boundary_condition("inlet", 2.0, 1);
                space.add_dirichlet_boundary_condition("inlet", 3.0, 2);

                space.add_dirichlet_boundary_condition("outlet", -1.0, 0);
                space.add_dirichlet_boundary_condition("outlet", -2.0, 1);
                space.add_dirichlet_boundary_condition("outlet", -3.0, 2);
            }

            VectorLaplaceOperator<3, Scalar_t> lapl{1.0};

            assemble_and_solve("vector_poisson", space, lapl);
        }

        void elasticity_problem() {
            static const int Dim = 3;
            auto params = input_params(get_more_complex_mesh_path(), Dim);

            FunctionSpace space;
            space.read(params);

            if (use_cube) {
                add_cube_bc(space, 3);
            } else {
                space.add_dirichlet_boundary_condition("inlet", -0.005, 0);
                space.add_dirichlet_boundary_condition("inlet", 0.0, 1);
                space.add_dirichlet_boundary_condition("inlet", 0.001, 2);

                space.add_dirichlet_boundary_condition("outlet", 0.005, 0);
                space.add_dirichlet_boundary_condition("outlet", 0.0, 1);
                space.add_dirichlet_boundary_condition("outlet", -0.001, 2);
            }

            LinearElasticity<Dim, Scalar_t> linear_elasticity{1.0, 1.0};
            assemble_and_solve("linear_elasticity", space, linear_elasticity);

            NeoHookean<Dim, Scalar_t> neohookean{1.0, 1.0};
            assemble_and_solve("neohookean", space, neohookean);
        }

        void poisson_problem_parallel_2D() {
            auto params = input_params(get_2D_mesh_path());

            FunctionSpace space;
            space.read(params);

            if (use_cube) {
                add_cube_bc(space, 1);
            } else {
                space.add_dirichlet_boundary_condition("surface_1", 1.0);
                space.add_dirichlet_boundary_condition("surface_3", -1.0);
            }

            LaplaceOperator<Scalar_t> lapl{1.0};

            std::stringstream ss;
            space.describe(ss);

            assemble_and_solve("poisson_problem_parallel_2D", space, lapl);
        }

        void poisson_problem_parallel_3D() {
            auto params = input_params(get_cube_path());

            FunctionSpace space;
            space.read(params);

            if (use_cube) {
                add_cube_bc(space, 1);
            } else {
                space.add_dirichlet_boundary_condition("body_top", 1.0);
                space.add_dirichlet_boundary_condition("body_bottom", -1.0);
            }

            LaplaceOperator<Scalar_t> lapl{1.0};

            std::stringstream ss;
            space.describe(ss);

            assemble_and_solve("poisson_problem_parallel_3D", space, lapl);
        }

        void elasticity_problem_parallel() {
            static const int Dim = 3;
            auto params = input_params(get_cube_path(), Dim);

            FunctionSpace space;
            space.read(params);

            if (use_cube) {
                add_cube_bc(space, 3);
            } else {
                space.add_dirichlet_boundary_condition("body_top", 0.0, 0);
                space.add_dirichlet_boundary_condition("body_top", -0.1, 1);
                space.add_dirichlet_boundary_condition("body_top", 0.0, 2);

                space.add_dirichlet_boundary_condition("body_bottom", 0.0, 0);
                space.add_dirichlet_boundary_condition("body_bottom", 0.1, 1);
                space.add_dirichlet_boundary_condition("body_bottom", 0.0, 2);
            }

            LinearElasticity<Dim, Scalar_t> linear_elasticity{1.0, 1.0};
            assemble_and_solve("elasticity_problem_parallel", space, linear_elasticity);
        }

        void shell_laplace_problem() {
            auto params = input_params(get_shell_mesh_path());

            FunctionSpace space;
            space.read(params);

            if (use_cube) {
                add_cube_bc(space, 1);
            } else {
                space.add_dirichlet_boundary_condition("Top", 1.0);
                space.add_dirichlet_boundary_condition("Side", -1.0);
            }

            LaplaceOperator<Scalar_t> lapl{1.0};

            std::stringstream ss;
            space.describe(ss);

            assemble_and_solve("shell_laplace_problem", space, lapl);
        }

        void shell_integral() {
            if (!use_cube) {
                auto params = input_params(get_shell_mesh_path());

                FunctionSpace space;
                space.read(params);

                auto fe_ptr = std::make_shared<FE>();
                create_fe(space, *fe_ptr, 0);
            } else {
                // FIXME
            }

            // TODO integrate whole surface

            // fe_ptr->print_jacobian();
            // fe_ptr->print_jacobian_inverse();
            // fe_ptr->print_measure();
        }

        void boundary_integral() {
            if (!use_cube) {
                auto params = input_params(get_2D_mesh_path());

                FunctionSpace space;
                space.read(params);

                auto fe_ptr = std::make_shared<FE>();
                create_fe_on_boundary(space, *fe_ptr, 0);
            }

            // TODO integrate whole marked boundary

            // fe_ptr->print_jacobian();
            // fe_ptr->print_jacobian_inverse();
            // fe_ptr->print_measure();
        }

        void run() {
            if (mpi_world_size() <= 4) {
                UTOPIA_RUN_TEST(create_fe_test);
                UTOPIA_RUN_TEST(poisson_problem);
                UTOPIA_RUN_TEST(vector_poisson_problem);
                save_output = export_tensors = true;
                UTOPIA_RUN_TEST(elasticity_problem);
                save_output = export_tensors = false;

                // FIXME for libmesh -> sideset nomenclature
                UTOPIA_RUN_TEST(poisson_problem_parallel_2D);
                UTOPIA_RUN_TEST(poisson_problem_parallel_3D);
                UTOPIA_RUN_TEST(elasticity_problem_parallel);
                UTOPIA_RUN_TEST(shell_integral);
                UTOPIA_RUN_TEST(shell_laplace_problem);
                UTOPIA_RUN_TEST(boundary_integral);
            }
        }
    };  // namespace utopia

}  // namespace utopia

#endif  // UTOPIA_SPACE_AND_FE_TEST_HPP
