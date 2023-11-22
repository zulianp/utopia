#include "PhaseField.h"
#include "PhaseFieldCoupled.h"

#include <basix/finite-element.h>
#include <climits>
#include <cmath>
#include <dolfinx.h>
#include <dolfinx/common/log.h>
#include <dolfinx/fem/assembler.h>
#include <dolfinx/fem/petsc.h>
#include <dolfinx/io/XDMFFile.h>
#include <dolfinx/la/Vector.h>
#include <dolfinx/mesh/Mesh.h>
#include <dolfinx/mesh/cell_types.h>

#include "utopia.hpp"

#include "utopia_TwoFieldAlternateMinimization.hpp"
#include "utopia_TwoFieldSPIN.hpp"

#include "DolfinxFunction.hpp"

using namespace dolfinx;
using T = PetscScalar;

std::shared_ptr<DolfinxFunction> phase_field_displacement(
    const std::shared_ptr<fem::FunctionSpace<T>> &V,
    std::shared_ptr<dolfinx::fem::Function<T>> u,
    std::shared_ptr<dolfinx::fem::Function<T>> c,
    std::vector<std::shared_ptr<const fem::DirichletBC<T>>> boundary_conditions) {
    auto objective = std::make_shared<fem::Form<T>>(
        fem::create_form<T>(*form_PhaseField_disp_objective, {}, {{"u", u}, {"c", c}}, {}, {}, V->mesh()));

    auto gradient = std::make_shared<fem::Form<T>>(
        fem::create_form<T>(*form_PhaseField_disp_gradient, {V}, {{"u", u}, {"c", c}}, {}, {}));

    auto hessian = std::make_shared<fem::Form<T>>(
        fem::create_form<T>(*form_PhaseField_disp_hessian, {V, V}, {{"u", u}, {"c", c}}, {}, {}));

    return std::make_shared<DolfinxFunction>(V, u, objective, gradient, hessian, boundary_conditions);
}

std::shared_ptr<DolfinxFunction> phase_field_phase(
    const std::shared_ptr<fem::FunctionSpace<T>> &V,
    std::shared_ptr<dolfinx::fem::Function<T>> u,
    std::shared_ptr<dolfinx::fem::Function<T>> c,
    std::vector<std::shared_ptr<const fem::DirichletBC<T>>> boundary_conditions) {
    auto objective = std::make_shared<fem::Form<T>>(
        fem::create_form<T>(*form_PhaseField_phase_objective, {}, {{"c", c}, {"u", u}}, {}, {}, V->mesh()));

    auto gradient = std::make_shared<fem::Form<T>>(
        fem::create_form<T>(*form_PhaseField_phase_gradient, {V}, {{"c", c}, {"u", u}}, {}, {}));

    auto hessian = std::make_shared<fem::Form<T>>(
        fem::create_form<T>(*form_PhaseField_phase_hessian, {V, V}, {{"c", c}, {"u", u}}, {}, {}));

    return std::make_shared<DolfinxFunction>(V, c, objective, gradient, hessian, boundary_conditions);
}

std::shared_ptr<DolfinxFunction> phase_field_coupled(
    const std::shared_ptr<fem::FunctionSpace<T>> &V,
    std::shared_ptr<dolfinx::fem::Function<T>> x,
    std::vector<std::shared_ptr<const fem::DirichletBC<T>>> boundary_conditions) {
    auto objective = std::make_shared<fem::Form<T>>(
        fem::create_form<T>(*form_PhaseFieldCoupled_objective, {}, {{"x", x}}, {}, {}, V->mesh()));

    auto gradient =
        std::make_shared<fem::Form<T>>(fem::create_form<T>(*form_PhaseFieldCoupled_gradient, {V}, {{"x", x}}, {}, {}));

    auto hessian = std::make_shared<fem::Form<T>>(
        fem::create_form<T>(*form_PhaseFieldCoupled_hessian, {V, V}, {{"x", x}}, {}, {}));

    return std::make_shared<DolfinxFunction>(V, x, objective, gradient, hessian, boundary_conditions);
}

auto disp_BC(const std::shared_ptr<fem::FunctionSpace<T>> &V) {
    auto bdofs_left = fem::locate_dofs_geometrical(*V, [](auto x) -> std::vector<std::int8_t> {
        constexpr T eps = 1.0e-6;
        std::vector<std::int8_t> marker(x.extent(1), false);
        for (std::size_t p = 0; p < x.extent(1); ++p) {
            if (std::abs(x(0, p)) < eps) marker[p] = true;
        }
        return marker;
    });

    auto bdofs_right = fem::locate_dofs_geometrical(*V, [](auto x) -> std::vector<std::int8_t> {
        constexpr T eps = 1.0e-6;
        std::vector<std::int8_t> marker(x.extent(1), false);
        for (std::size_t p = 0; p < x.extent(1); ++p) {
            if (std::abs(x(0, p) - 1) < eps) marker[p] = true;
        }
        return marker;
    });

    return std::vector{std::make_shared<const fem::DirichletBC<T>>(std::vector<T>{-1e-5, 0, 0}, bdofs_left, V),
                       std::make_shared<const fem::DirichletBC<T>>(std::vector<T>{0, 0, 0}, bdofs_right, V)};
}

auto phase_BC(const std::shared_ptr<fem::FunctionSpace<T>> &C) {
    auto frac_locator = fem::locate_dofs_geometrical(*C, [](auto x) -> std::vector<std::int8_t> {
        std::vector<std::int8_t> marker(x.extent(1), false);
        for (std::size_t p = 0; p < x.extent(1); ++p) {
            marker[p] = 0.4 <= x(0, p) && x(0, p) <= 0.6 && x(1, p) <= 0.5;
        }
        return marker;
    });

    return std::vector{std::make_shared<const fem::DirichletBC<T>>(1., frac_locator, C)};
    // return std::vector{std::make_shared<const fem::DirichletBC<T>>(0., frac_locator, C)};
}

// auto coupled_BC(const std::shared_ptr<fem::FunctionSpace> &V)
// {
//     auto bdofs_left = fem::locate_dofs_geometrical(
//         {*V}, [](auto &&x) -> xt::xtensor<bool, 1> { return xt::isclose(xt::row(x, 0), 0.0); });

//     auto bdofs_right = fem::locate_dofs_geometrical(
//         {*V}, [](auto &&x) -> xt::xtensor<bool, 1> { return xt::isclose(xt::row(x, 0), 1.0); });

//     auto frac_locator = fem::locate_dofs_geometrical({*V}, [](auto &&p) -> xt::xtensor<bool, 1> {
//         auto x = xt::row(p, 0);
//         auto y = xt::row(p, 1);
//         auto z = xt::row(p, 2);
//         return 0.49 <= x <= 0.51 && y <= 0.5;
//     });

//     return std::vector{
//         std::make_shared<const fem::DirichletBC<T>>(xt::xarray<T>{0, 0, 0, 0}, bdofs_left, V),
//         std::make_shared<const fem::DirichletBC<T>>(xt::xarray<T>{1e-5, 0, 0, 0}, bdofs_right, V),
//         std::make_shared<const fem::DirichletBC<T>>(xt::xarray<T>{1e-5, 0, 0, 1}, frac_locator, V)};
// }

int main(int argc, char *argv[]) {
    dolfinx::init_logging(argc, argv);
    utopia::Utopia::Init(argc, argv);

    // Set the logging thread name to show the process rank
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    std::string thread_name = "RANK " + std::to_string(mpi_rank);
    loguru::set_thread_name(thread_name.c_str());

    {
        static int dims = 3;
        auto mesh =
            std::make_shared<mesh::Mesh<T>>(mesh::create_box(MPI_COMM_WORLD,
                                                          {{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}}},
                                                          {10, 10, 11},
                                                          mesh::CellType::tetrahedron,
                                                          mesh::create_cell_partitioner(mesh::GhostMode::none)));

        auto V = std::make_shared<fem::FunctionSpace<T>>(
            fem::create_functionspace(functionspace_form_PhaseField_disp_gradient, "u", mesh));

        auto C = std::make_shared<fem::FunctionSpace<T>>(
            fem::create_functionspace(functionspace_form_PhaseField_phase_gradient, "c", mesh));

        auto X = std::make_shared<fem::FunctionSpace<T>>(
            fem::create_functionspace(functionspace_form_PhaseFieldCoupled_gradient, "x", mesh));

        auto u = std::make_shared<fem::Function<T>>(V);
        auto c = std::make_shared<fem::Function<T>>(C);
        auto x = std::make_shared<fem::Function<T>>(X);

        ////////////////////////////////////////////////////////////////
        // Create Dirichlet boundary conditions
        ////////////////////////////////////////////////////////////////

        auto disp_bcs = disp_BC(V);
        auto phase_bcs = phase_BC(C);

        auto X_disp = X->sub({0})->collapse();
        auto X_phase = X->sub({1})->collapse();

        auto coupled_bcs = disp_BC(utopia::make_ref(X_disp.first));
        auto X_phase_bcs = phase_BC(utopia::make_ref(X_phase.first));
        coupled_bcs.insert(coupled_bcs.end(), X_phase_bcs.begin(), X_phase_bcs.end());

        // auto coupled_bcs = coupled_BC(X);

        ////////////////////////////////////////////////////////////////

        bool verbose = true;

        // Split-fields
        auto f1 = phase_field_displacement(V, u, c, disp_bcs);
        auto f2 = phase_field_phase(C, u, c, phase_bcs);

        // Global optimization problem (coupled fields)
        auto c12 = phase_field_coupled(X, x, coupled_bcs);

        auto nls1 = std::make_shared<utopia::Newton<utopia::PetscMatrix>>(
            std::make_shared<utopia::BiCGStab<utopia::PetscMatrix, utopia::PetscVector, utopia::HOMEMADE>>());

        auto nls2 = std::make_shared<utopia::Newton<utopia::PetscMatrix>>(
            std::make_shared<utopia::BiCGStab<utopia::PetscMatrix, utopia::PetscVector, utopia::HOMEMADE>>());

        // nls1->linear_solver()->verbose(true);
        // nls2->linear_solver()->verbose(true);

        // nls1->verbose(verbose);
        // nls2->verbose(verbose);

        auto f1_to_c12 = [](const utopia::PetscVector &in, utopia::PetscVector &out) {
            // utopia::out() << "f1_to_c12: " << in.size() << " -> " << out.size() << "\n";

            using namespace utopia;

            auto in_view = local_view_device(in);
            auto out_view = local_view_device(out);

            parallel_for(
                local_range_device(in), UTOPIA_LAMBDA(SizeType i) {
                    SizeType ii = i / dims;
                    int rmd = i - ii * dims;
                    auto v = in_view.get(i);

                    auto j = ii * (dims + 1) + rmd;
                    out_view.set(j, v);
                });
        };

        auto f2_to_c12 = [](const utopia::PetscVector &in, utopia::PetscVector &out) {
            // utopia::out() << "f2_to_c12: " << in.size() << " -> " << out.size() << "\n";

            using namespace utopia;

            auto in_view = local_view_device(in);
            auto out_view = local_view_device(out);

            parallel_for(
                local_range_device(in), UTOPIA_LAMBDA(SizeType i) {
                    auto v = in_view.get(i);
                    out_view.set(i * (dims + 1) + dims, v);
                });
        };

        auto c12_to_f1 = [](const utopia::PetscVector &in, utopia::PetscVector &out) {
            // utopia::out() << "c12_to_f1: " << in.size() << " -> " << out.size() << "\n";

            using namespace utopia;

            auto in_view = local_view_device(in);
            auto out_view = local_view_device(out);

            parallel_for(
                local_range_device(out), UTOPIA_LAMBDA(SizeType i) {
                    SizeType ii = i / dims;
                    int rmd = i - ii * dims;
                    auto v = in_view.get(ii * (dims + 1) + rmd);
                    out_view.set(i, v);
                });
        };

        auto c12_to_f2 = [](const utopia::PetscVector &in, utopia::PetscVector &out) {
            // utopia::out() << "c12_to_f2: " << in.size() << " -> " << out.size() << "\n";

            using namespace utopia;

            auto in_view = local_view_device(in);
            auto out_view = local_view_device(out);

            parallel_for(
                local_range_device(out), UTOPIA_LAMBDA(SizeType i) {
                    auto v = in_view.get(i * (dims + 1) + dims);
                    out_view.set(i, v);
                });
        };

        utopia::PetscVector solution;
        c12->create_vector(solution);
        solution.set(1);

        if (true) {
            // auto ls = std::make_shared<utopia::BiCGStab<utopia::PetscMatrix, utopia::PetscVector,
            // utopia::HOMEMADE>>(); ls->verbose(true);

            auto ls = std::make_shared<utopia::Factorization<utopia::PetscMatrix, utopia::PetscVector>>();

            utopia::Newton<utopia::PetscMatrix> newton(ls);

            auto params = utopia::param_list(   //
                utopia::param("damping", 0.5),  //
                utopia::param("max_it", 0),     //
                utopia::param("verbose", true)  //
            );

            newton.read(params);
            // newton.solve(*c12, solution);

        } else if   //
            (true)  //
        // (false)  //
        {
            utopia::TwoFieldAlternateMinimization<utopia::PetscMatrix> tfa(nls1, nls2);
            tfa.set_field_functions(f1, f2);
            tfa.set_transfers(f1_to_c12, f2_to_c12, c12_to_f1, c12_to_f2);
            tfa.verbose(verbose);
            tfa.max_it(400);
            // tfa.verbosity_level(utopia::VERBOSITY_LEVEL_DEBUG);
            // tfa.verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);

            tfa.field1_diff_tol(1e-14);
            tfa.field2_diff_tol(1e-14);
            tfa.solve(*c12, solution);
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
            tfs.solve(*c12, solution);

            /////////////////////////////////////////////////////////
            // Save solution in VTK format
        }

        {
            // Make sure that we plot the global solution
            utopia::PetscVector c12_temp, f1_temp, f2_temp;

            f1->create_vector(f1_temp);
            f2->create_vector(f2_temp);
            c12->create_vector(c12_temp);

            c12_to_f1(c12_temp, f1_temp);
            c12_to_f2(c12_temp, f2_temp);
        }

        io::VTKFile file_u(mesh->comm(), "u.pvd", "w");
        file_u.write<T>({*f1->u()}, 0.0);

        io::VTKFile file_c(mesh->comm(), "c.pvd", "w");
        file_c.write<T>({*f2->u()}, 0.0);
    }

    // PetscFinalize();
    return utopia::Utopia::Finalize();
}
