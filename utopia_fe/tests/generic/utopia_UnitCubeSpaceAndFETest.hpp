#ifndef UTOPIA_UNIT_CUBE_SPACE_AND_FE_TEST_HPP
#define UTOPIA_UNIT_CUBE_SPACE_AND_FE_TEST_HPP

#include "utopia_Testing.hpp"
#include "utopia_fe_base.hpp"
#include "utopia_ui.hpp"

#include "utopia_fe_Core.hpp"

// FIXME
#include "utopia_kokkos_LaplaceOperator.hpp"
#include "utopia_kokkos_LinearElasticity.hpp"
#include "utopia_kokkos_NeoHookean.hpp"
#include "utopia_kokkos_VectorLaplaceOperator.hpp"

namespace utopia {

    template <class FunctionSpace, class FE>
    class UnitCubeSpaceAndFETest {
    public:
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using SideSet_t = typename Traits<FunctionSpace>::SideSet;
        using Field_t = utopia::Field<FunctionSpace>;

        bool save_output{false};
        bool export_tensors{false};
        bool verbose{false};
        std::string output_format{"e"};
        Scalar_t rtol{1e-6};

        InputParameters cube_space_param(const int n_var) const {
            const int nx = 10;
            const int ny = 10;
            const int nz = 10;

            return param_list(
                param("n_var", n_var),
                param("mesh", param_list(param("type", "cube"), param("nx", nx), param("ny", ny), param("nz", nz))));
        }

        inline void add_cube_bc(FunctionSpace &space, const int n_var = 1) const {
            for (int i = 0; i < n_var; ++i) {
                space.add_dirichlet_boundary_condition(SideSet_t::Cube::left(), -(i == 0), i);
                space.add_dirichlet_boundary_condition(SideSet_t::Cube::right(), (i == 0), i);
            }
        }

        template <class Assembler>
        void assemble_and_solve(const std::string &name, FunctionSpace &space, Assembler &assembler) {
            using FEField_t = typename Traits<FE>::Field;

            Vector_t rhs;
            space.create_vector(rhs);

            Field_t field_x;
            space.create_field(field_x);

            auto &x = field_x.data();
            x.set(0.0);

            FEField_t fe_field(assembler.fe_ptr());
            convert_field(field_x, fe_field);

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
            solver.rtol(rtol);
            solver.verbose(verbose);

            utopia_test_assert(solver.solve(mat, rhs, x));

            ////////////////////////////////////////////////////////////

            if (save_output) {
                space.write(name + "." + output_format, x);
                space.write("rhs." + output_format, rhs);
            }

            if (export_tensors) {
                write(name + "_x.m", x);
            }
        }

        virtual void run() {
            UTOPIA_RUN_TEST(unit_cube_create_fe);

            save_output = export_tensors = true;
            UTOPIA_RUN_TEST(unit_cube_poisson_problem);
            save_output = export_tensors = false;

            UTOPIA_RUN_TEST(unit_cube_vector_poisson_problem);
            UTOPIA_RUN_TEST(unit_cube_elasticity_problem);
            UTOPIA_RUN_TEST(unit_cube_poisson_problem_parallel_2D);
            UTOPIA_RUN_TEST(unit_cube_poisson_problem_parallel_3D);
            UTOPIA_RUN_TEST(unit_cube_elasticity_problem_parallel);
        }

    private:
        void unit_cube_create_fe() {
            auto params = cube_space_param(1);
            FunctionSpace space;
            space.read(params);

            auto fe_ptr = std::make_shared<FE>();
            create_fe(space, *fe_ptr, 0);
        }

        void unit_cube_poisson_problem() {
            auto params = cube_space_param(1);
            FunctionSpace space;
            space.read(params);

            add_cube_bc(space, 1);

            auto fe_ptr = std::make_shared<FE>();
            create_fe(space, *fe_ptr, 2);
            utopia::kokkos::LaplaceOperator<FE> lapl(fe_ptr, {1.0});

            assemble_and_solve("poisson", space, lapl);
        }

        void unit_cube_vector_poisson_problem() {
            auto params = cube_space_param(3);

            FunctionSpace space;
            space.read(params);

            add_cube_bc(space, 3);

            auto fe_ptr = std::make_shared<FE>();
            create_fe(space, *fe_ptr, 2);
            utopia::kokkos::VectorLaplaceOperator<FE, 3> lapl(fe_ptr, {1.0});

            assemble_and_solve("vector_poisson", space, lapl);
        }

        void unit_cube_elasticity_problem() {
            static const int Dim = 3;
            auto params = cube_space_param(Dim);

            FunctionSpace space;
            space.read(params);

            add_cube_bc(space, 3);

            auto fe_ptr = std::make_shared<FE>();
            create_fe(space, *fe_ptr, 2);

            utopia::kokkos::LinearElasticity<FE, Dim> linear_elasticity(fe_ptr, {1.0, 1.0});
            assemble_and_solve("linear_elasticity", space, linear_elasticity);

            utopia::kokkos::NeoHookean<FE> neohookean(fe_ptr, {1.0, 1.0});
            assemble_and_solve("neohookean", space, neohookean);
        }

        void unit_cube_poisson_problem_parallel_2D() {
            auto params = cube_space_param(1);

            FunctionSpace space;
            space.read(params);

            add_cube_bc(space, 1);

            auto fe_ptr = std::make_shared<FE>();
            create_fe(space, *fe_ptr, 2);
            utopia::kokkos::LaplaceOperator<FE> lapl(fe_ptr, {1.0});

            std::stringstream ss;
            space.describe(ss);

            assemble_and_solve("unit_cube_poisson_problem_parallel_2D", space, lapl);
        }

        void unit_cube_poisson_problem_parallel_3D() {
            auto params = cube_space_param(1);

            FunctionSpace space;
            space.read(params);

            add_cube_bc(space, 1);

            auto fe_ptr = std::make_shared<FE>();
            create_fe(space, *fe_ptr, 2);
            utopia::kokkos::LaplaceOperator<FE> lapl(fe_ptr, {1.0});

            std::stringstream ss;
            space.describe(ss);

            assemble_and_solve("unit_cube_poisson_problem_parallel_3D", space, lapl);
        }

        void unit_cube_elasticity_problem_parallel() {
            static const int Dim = 3;
            auto params = cube_space_param(Dim);

            FunctionSpace space;
            space.read(params);

            add_cube_bc(space, 3);

            auto fe_ptr = std::make_shared<FE>();
            create_fe(space, *fe_ptr, 2);

            utopia::kokkos::LinearElasticity<FE, Dim> linear_elasticity(fe_ptr, {1.0, 1.0});

            assemble_and_solve("unit_cube_elasticity_problem_parallel", space, linear_elasticity);
        }
    };
}  // namespace utopia

#endif  // UTOPIA_UNIT_CUBE_SPACE_AND_FE_TEST_HPP
