#ifndef UTOPIA_SPACE_AND_FE_TEST_HPP
#define UTOPIA_SPACE_AND_FE_TEST_HPP

#include "utopia_Testing.hpp"
#include "utopia_fe_base.hpp"
#include "utopia_ui.hpp"

#include "utopia_fe_Core.hpp"

// FIXME
#include "utopia_kokkos_LaplaceOperator.hpp"
#include "utopia_kokkos_LinearElasticity.hpp"
#include "utopia_kokkos_NeoHookean.hpp"
#include "utopia_kokkos_VectorLaplaceOperator.hpp"

#include "utopia_UnitCubeSpaceAndFETest.hpp"
namespace utopia {

    template <class FunctionSpace, class FE>
    class SpaceAndFETest : public UnitCubeSpaceAndFETest<FunctionSpace, FE> {
    public:
        using Super = utopia::UnitCubeSpaceAndFETest<FunctionSpace, FE>;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using SideSet_t = typename Traits<FunctionSpace>::SideSet;
        using Field_t = utopia::Field<FunctionSpace>;

        bool test_data_unavailable{true};

        static std::string get_shell_mesh_path() { return "../data/fe_problem_solve/shell.e"; }

        static std::string get_cube_path() {
            std::string dir = "../data/knf/cube_vs_cube";
            if (mpi_world_size() > 1) {
                dir += "/" + std::to_string(mpi_world_size());
            }

            return dir + "/body.e";
        }

        InputParameters input_params(const std::string &path, const int n_var = 1) const {
            if (test_data_unavailable) {
                return this->cube_space_param(n_var);
            } else {
                return param_list(param("n_var", 3),
                                  param("mesh", param_list(param("type", "file"), param("path", path))));
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

        void poisson_problem() {
            auto params = input_params(get_more_complex_mesh_path());
            FunctionSpace space;
            space.read(params);

            if (test_data_unavailable) {
                this->add_cube_bc(space, 1);
            } else {
                space.add_dirichlet_boundary_condition("inlet", 1.0);
                space.add_dirichlet_boundary_condition("outlet", -1.0);
            }

            auto fe_ptr = std::make_shared<FE>();
            create_fe(space, *fe_ptr, 2);
            utopia::kokkos::LaplaceOperator<FunctionSpace, FE> lapl(fe_ptr, {1.0});

            Super::assemble_and_solve("poisson", space, lapl);
        }

        void vector_poisson_problem() {
            auto params = input_params(get_more_complex_mesh_path(), 3);

            FunctionSpace space;
            space.read(params);

            if (test_data_unavailable) {
                this->add_cube_bc(space, 3);
            } else {
                space.add_dirichlet_boundary_condition("inlet", 1.0, 0);
                space.add_dirichlet_boundary_condition("inlet", 2.0, 1);
                space.add_dirichlet_boundary_condition("inlet", 3.0, 2);

                space.add_dirichlet_boundary_condition("outlet", -1.0, 0);
                space.add_dirichlet_boundary_condition("outlet", -2.0, 1);
                space.add_dirichlet_boundary_condition("outlet", -3.0, 2);
            }

            auto fe_ptr = std::make_shared<FE>();
            create_fe(space, *fe_ptr, 2);
            utopia::kokkos::VectorLaplaceOperator<FunctionSpace, FE, 3> lapl(fe_ptr, {1.0});

            Super::assemble_and_solve("vector_poisson", space, lapl);
        }

        void elasticity_problem() {
            static const int Dim = 3;
            auto params = input_params(get_more_complex_mesh_path(), Dim);

            FunctionSpace space;
            space.read(params);

            if (test_data_unavailable) {
                this->add_cube_bc(space, 3);
            } else {
                space.add_dirichlet_boundary_condition("inlet", -0.005, 0);
                space.add_dirichlet_boundary_condition("inlet", 0.0, 1);
                space.add_dirichlet_boundary_condition("inlet", 0.001, 2);

                space.add_dirichlet_boundary_condition("outlet", 0.005, 0);
                space.add_dirichlet_boundary_condition("outlet", 0.0, 1);
                space.add_dirichlet_boundary_condition("outlet", -0.001, 2);
            }

            auto fe_ptr = std::make_shared<FE>();
            create_fe(space, *fe_ptr, 2);

            utopia::kokkos::LinearElasticity<FunctionSpace, FE, Dim> linear_elasticity(fe_ptr, {1.0, 1.0});
            Super::assemble_and_solve("linear_elasticity", space, linear_elasticity);

            utopia::kokkos::NeoHookean<FunctionSpace, FE> neohookean(fe_ptr, {1.0, 1.0});
            Super::assemble_and_solve("neohookean", space, neohookean);
        }

        void poisson_problem_parallel_2D() {
            auto params = input_params(get_2D_mesh_path());

            FunctionSpace space;
            space.read(params);

            if (test_data_unavailable) {
                this->add_cube_bc(space, 1);
            } else {
                space.add_dirichlet_boundary_condition("surface_1", 1.0);
                space.add_dirichlet_boundary_condition("surface_3", -1.0);
            }

            auto fe_ptr = std::make_shared<FE>();
            create_fe(space, *fe_ptr, 2);
            utopia::kokkos::LaplaceOperator<FunctionSpace, FE> lapl(fe_ptr, {1.0});

            std::stringstream ss;
            space.describe(ss);

            Super::assemble_and_solve("poisson_problem_parallel_2D", space, lapl);
        }

        void poisson_problem_parallel_3D() {
            auto params = input_params(get_cube_path());

            FunctionSpace space;
            space.read(params);

            if (test_data_unavailable) {
                this->add_cube_bc(space, 1);
            } else {
                space.add_dirichlet_boundary_condition("body_top", 1.0);
                space.add_dirichlet_boundary_condition("body_bottom", -1.0);
            }

            auto fe_ptr = std::make_shared<FE>();
            create_fe(space, *fe_ptr, 2);
            utopia::kokkos::LaplaceOperator<FunctionSpace, FE> lapl(fe_ptr, {1.0});

            std::stringstream ss;
            space.describe(ss);

            Super::assemble_and_solve("poisson_problem_parallel_3D", space, lapl);
        }

        void elasticity_problem_parallel() {
            static const int Dim = 3;
            auto params = input_params(get_cube_path(), Dim);

            FunctionSpace space;
            space.read(params);

            if (test_data_unavailable) {
                this->add_cube_bc(space, 3);
            } else {
                space.add_dirichlet_boundary_condition("body_top", 0.0, 0);
                space.add_dirichlet_boundary_condition("body_top", -0.1, 1);
                space.add_dirichlet_boundary_condition("body_top", 0.0, 2);

                space.add_dirichlet_boundary_condition("body_bottom", 0.0, 0);
                space.add_dirichlet_boundary_condition("body_bottom", 0.1, 1);
                space.add_dirichlet_boundary_condition("body_bottom", 0.0, 2);
            }

            auto fe_ptr = std::make_shared<FE>();
            create_fe(space, *fe_ptr, 2);

            utopia::kokkos::LinearElasticity<FunctionSpace, FE, Dim> linear_elasticity(fe_ptr, {1.0, 1.0});

            Super::assemble_and_solve("elasticity_problem_parallel", space, linear_elasticity);
        }

        void shell_laplace_problem() {
            auto params = input_params(get_shell_mesh_path());

            FunctionSpace space;
            space.read(params);

            if (test_data_unavailable) {
                this->add_cube_bc(space, 1);
            } else {
                space.add_dirichlet_boundary_condition("Top", 1.0);
                space.add_dirichlet_boundary_condition("Side", -1.0);
            }

            auto fe_ptr = std::make_shared<FE>();
            create_fe(space, *fe_ptr, 2);
            utopia::kokkos::LaplaceOperator<FunctionSpace, FE> lapl(fe_ptr, {1.0});

            std::stringstream ss;
            space.describe(ss);

            Super::assemble_and_solve("shell_laplace_problem", space, lapl);
        }

        void shell_integral() {
            if (!test_data_unavailable) {
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
            if (!test_data_unavailable) {
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
            Super::run();

            if (mpi_world_size() <= 4) {
                UTOPIA_RUN_TEST(poisson_problem);
                UTOPIA_RUN_TEST(vector_poisson_problem);
                this->save_output = this->export_tensors = true;
                UTOPIA_RUN_TEST(elasticity_problem);
                this->save_output = this->export_tensors = false;

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
