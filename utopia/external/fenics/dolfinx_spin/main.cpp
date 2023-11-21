#include <basix/finite-element.h>
#include <dolfinx.h>
#include <dolfinx/common/log.h>
#include <dolfinx/fem/assembler.h>
#include <dolfinx/fem/petsc.h>
#include <dolfinx/io/XDMFFile.h>
#include <dolfinx/la/Vector.h>
#include <dolfinx/mesh/Mesh.h>
#include <dolfinx/mesh/cell_types.h>
#include <dolfinx/nls/NewtonSolver.h>
#include <cmath>
#include <xtensor/xarray.hpp>
#include <xtensor/xview.hpp>

#include "PhaseField.h"
#include "PhaseFieldCoupled.h"

#include "utopia.hpp"

#include "utopia_TwoFieldAlternateMinimization.hpp"
#include "utopia_TwoFieldSPIN.hpp"

#include "DolfinxFunction.hpp"

using namespace dolfinx;
using T = PetscScalar;

std::shared_ptr<DolfinxFunction> phase_field_displacement(
    const std::shared_ptr<fem::FunctionSpace> &V,
    std::vector<std::shared_ptr<const fem::DirichletBC<T>>> boundary_conditions) {
    auto u = std::make_shared<fem::Function<T>>(V);

    auto objective = std::make_shared<fem::Form<T>>(
        fem::create_form<T>(*form_PhaseField_disp_objective, {}, {{"u", u}}, {}, {}, V->mesh()));

    auto gradient = std::make_shared<fem::Form<T>>(
        fem::create_form<T>(*form_PhaseField_disp_gradient, {V}, {{"u", u}}, {}, {}));

    auto hessian = std::make_shared<fem::Form<T>>(
        fem::create_form<T>(*form_PhaseField_disp_hessian, {V, V}, {{"u", u}}, {}, {}));

    return std::make_shared<DolfinxFunction>(V, u, objective, gradient, hessian, boundary_conditions);
}

std::shared_ptr<DolfinxFunction> phase_field_phase(
    const std::shared_ptr<fem::FunctionSpace> &V,
    std::vector<std::shared_ptr<const fem::DirichletBC<T>>> boundary_conditions) {
    auto u = std::make_shared<fem::Function<T>>(V);

    auto objective = std::make_shared<fem::Form<T>>(
        fem::create_form<T>(*form_PhaseField_phase_objective, {}, {{"u", u}}, {}, {}, V->mesh()));

    auto gradient = std::make_shared<fem::Form<T>>(
        fem::create_form<T>(*form_PhaseField_phase_gradient, {V}, {{"u", u}}, {}, {}));

    auto hessian = std::make_shared<fem::Form<T>>(
        fem::create_form<T>(*form_PhaseField_phase_hessian, {V, V}, {{"u", u}}, {}, {}));

    return std::make_shared<DolfinxFunction>(V, u, objective, gradient, hessian, boundary_conditions);
}

std::shared_ptr<DolfinxFunction> phase_field_coupled(
    const std::shared_ptr<fem::FunctionSpace> &V,
    std::vector<std::shared_ptr<const fem::DirichletBC<T>>> boundary_conditions) {
    auto u = std::make_shared<fem::Function<T>>(V);

    auto objective = std::make_shared<fem::Form<T>>(
        fem::create_form<T>(*form_PhaseFieldCoupled_objective, {}, {{"u", u}}, {}, {}, V->mesh()));

    auto gradient = std::make_shared<fem::Form<T>>(
        fem::create_form<T>(*form_PhaseFieldCoupled_gradient, {V}, {{"u", u}}, {}, {}));

    auto hessian = std::make_shared<fem::Form<T>>(
        fem::create_form<T>(*form_PhaseFieldCoupled_hessian, {V, V}, {{"u", u}}, {}, {}));

    return std::make_shared<DolfinxFunction>(V, u, objective, gradient, hessian, boundary_conditions);
}

int main(int argc, char *argv[]) {
    dolfinx::init_logging(argc, argv);
    utopia::Utopia::Init(argc, argv);

    // Set the logging thread name to show the process rank
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    std::string thread_name = "RANK " + std::to_string(mpi_rank);
    loguru::set_thread_name(thread_name.c_str());

    {
        // Inside the ``main`` function, we begin by defining a tetrahedral mesh
        // of the domain and the function space on this mesh. Here, we choose to
        // create a unit cube mesh with 25 ( = 24 + 1) vertices in one direction
        // and 17 ( = 16 + 1) vertices in the other two directions. With this
        // mesh, we initialize the (finite element) function space defined by the
        // generated code.
        //
        // .. code-block:: cpp

        // Create mesh and define function space
        auto mesh = std::make_shared<mesh::Mesh>(mesh::create_box(MPI_COMM_WORLD,
                                                                  {{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}}},
                                                                  {10, 10, 10},
                                                                  mesh::CellType::tetrahedron,
                                                                  mesh::GhostMode::none));

        auto V = std::make_shared<fem::FunctionSpace>(
            fem::create_functionspace(functionspace_form_PhaseField_disp_gradient, "u", mesh));

        auto u_rotation = std::make_shared<fem::Function<T>>(V);
        u_rotation->interpolate([](auto &&x) {
            constexpr double scale = 0.1;

            // Center of rotation
            constexpr double x1_c = 0.5;
            constexpr double x2_c = 0.5;

            // Large angle of rotation (60 degrees)
            constexpr double theta = 1.04719755;

            // New coordinates
            auto x1 = xt::row(x, 1);
            auto x2 = xt::row(x, 2);
            xt::xarray<double> values = xt::zeros_like(x);
            xt::row(values, 1) = scale * (x1_c + (x1 - x1_c) * std::cos(theta) - (x2 - x2_c) * std::sin(theta) - x1);
            xt::row(values, 2) = scale * (x2_c + (x1 - x1_c) * std::sin(theta) - (x2 - x2_c) * std::cos(theta) - x2);
            return values;
        });

        // Create Dirichlet boundary conditions
        auto bdofs_left = fem::locate_dofs_geometrical(
            {*V}, [](auto &&x) -> xt::xtensor<bool, 1> { return xt::isclose(xt::row(x, 0), 0.0); });

        auto bdofs_right = fem::locate_dofs_geometrical(
            {*V}, [](auto &&x) -> xt::xtensor<bool, 1> { return xt::isclose(xt::row(x, 0), 1.0); });

        auto bcs = std::vector{std::make_shared<const fem::DirichletBC<T>>(xt::xarray<T>{0, 0, 0}, bdofs_left, V),
                               std::make_shared<const fem::DirichletBC<T>>(u_rotation, bdofs_right)};

        /////////////////////////////////////////////////////////

        bool verbose = true;

        // Split-fields
        auto f1 = phase_field_displacement(V, bcs);
        auto f2 = phase_field_phase(V, bcs);

        // Global optimization problem (coupled fields)
        auto c12 = phase_field_coupled(V, bcs);

        auto nls1 = std::make_shared<utopia::Newton<utopia::PetscMatrix>>(
            std::make_shared<utopia::BiCGStab<utopia::PetscMatrix, utopia::PetscVector, utopia::HOMEMADE>>());

        auto nls2 = std::make_shared<utopia::Newton<utopia::PetscMatrix>>(
            std::make_shared<utopia::BiCGStab<utopia::PetscMatrix, utopia::PetscVector, utopia::HOMEMADE>>());

        nls1->linear_solver()->verbose(true);
        nls2->linear_solver()->verbose(true);

        nls1->verbose(verbose);
        nls2->verbose(verbose);

        auto f1_to_c12 = [](const utopia::PetscVector &in, utopia::PetscVector &out) { out = in; };
        auto f2_to_c12 = f1_to_c12;
        auto c12_to_f1 = f1_to_c12;
        auto c12_to_f2 = f1_to_c12;

        utopia::PetscVector x;
        c12->create_vector(x);

        if (0) {
            utopia::TwoFieldAlternateMinimization<utopia::PetscMatrix> tfa(nls1, nls2);
            tfa.set_field_functions(f1, f2);
            tfa.set_transfers(f1_to_c12, f2_to_c12, c12_to_f1, c12_to_f2);
            tfa.verbose(verbose);
            tfa.verbosity_level(utopia::VERBOSITY_LEVEL_DEBUG);
            tfa.solve(*c12, x);
        } else {
            auto lsc12 = std::make_shared<utopia::KSP_MF<utopia::PetscMatrix, utopia::PetscVector>>();
            lsc12->ksp_type("gmres");
            lsc12->pc_type("none");

            lsc12->verbose(verbose);
            nls1->verbose(verbose);
            nls2->verbose(verbose);

            utopia::TwoFieldSPIN<utopia::PetscMatrix> tfs(lsc12, nls1, nls2);
            tfs.set_field_functions(f1, f2);
            tfs.set_transfers(f1_to_c12, f2_to_c12, c12_to_f1, c12_to_f2);
            // tfs.additive_precond(false);
            tfs.additive_precond(true);
            tfs.verbosity_level(utopia::VERBOSITY_LEVEL_DEBUG);
            tfs.solve(*c12, x);
        }

        /////////////////////////////////////////////////////////
        // Save solution in VTK format
        io::VTKFile file_u(mesh->comm(), "u.pvd", "w");
        file_u.write<T>({*c12->u()}, 0.0);
    }

    // PetscFinalize();
    return utopia::Utopia::Finalize();
}
