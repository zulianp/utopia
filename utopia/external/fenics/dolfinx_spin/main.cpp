#include "PhaseField.h"
#include "PhaseFieldCoupled.h"

#include <basix/finite-element.h>
#include <dolfinx.h>
#include <dolfinx/common/log.h>
#include <dolfinx/fem/assembler.h>
#include <dolfinx/fem/petsc.h>
#include <dolfinx/io/XDMFFile.h>
#include <dolfinx/la/Vector.h>
#include <dolfinx/mesh/Mesh.h>
#include <dolfinx/mesh/cell_types.h>
#include <climits>
#include <cmath>

#include "utopia.hpp"

#include "utopia_PseudoTimeStepper.hpp"
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

auto create_frac_locator(const fem::FunctionSpace<T> &C) {
    return fem::locate_dofs_geometrical(C, [](auto x) -> std::vector<std::int8_t> {
        std::vector<std::int8_t> marker(x.extent(1), false);
        for (std::size_t p = 0; p < x.extent(1); ++p) {
            marker[p] = (x(0, p) >= 0.4) && (x(0, p) <= 0.6) && (x(1, p) <= 0.5);
        }
        return marker;
    });
}

auto create_left_boundary_locator(const fem::FunctionSpace<T> &V) {
    return fem::locate_dofs_geometrical(V, [](auto x) -> std::vector<std::int8_t> {
        constexpr T eps = 1.0e-6;
        std::vector<std::int8_t> marker(x.extent(1), false);
        for (std::size_t p = 0; p < x.extent(1); ++p) {
            if (std::abs(x(0, p)) < eps) marker[p] = true;
        }

        return marker;
    });
}

auto create_right_boundary_locator(const fem::FunctionSpace<T> &V) {
    return fem::locate_dofs_geometrical(V, [](auto x) -> std::vector<std::int8_t> {
        constexpr T eps = 1.0e-6;
        std::vector<std::int8_t> marker(x.extent(1), false);
        for (std::size_t p = 0; p < x.extent(1); ++p) {
            if (std::abs(x(0, p) - 1) < eps) marker[p] = true;
        }

        return marker;
    });
}

auto disp_BC(const std::shared_ptr<fem::FunctionSpace<T>> &V, const T disp_x) {
    auto bdofs_left = create_left_boundary_locator(*V);
    auto bdofs_right = create_right_boundary_locator(*V);

    return std::vector{std::make_shared<const fem::DirichletBC<T>>(std::vector<T>{0, 0, 0}, bdofs_left, V),
                       std::make_shared<const fem::DirichletBC<T>>(std::vector<T>{disp_x, 0, 0}, bdofs_right, V)};
}

auto phase_BC(const std::shared_ptr<fem::FunctionSpace<T>> &C) {
    auto frac_locator = create_frac_locator(*C);
    return std::vector{std::make_shared<const fem::DirichletBC<T>>(1., frac_locator, C)};
}

auto coupled_BC(const std::shared_ptr<fem::FunctionSpace<T>> &X, const T disp_x) {
    std::vector<std::shared_ptr<const fem::DirichletBC<T>>> bcs;

    static const int dim = 3;
    T imposed_disp[3] = {disp_x, 0, 0};

    // Displacement BCs
    for (int d = 0; d < dim; d++) {
        auto X_sub = X->sub({0})->sub({d});
        auto VxDofs = X_sub->collapse();
        auto &V = VxDofs.first;

        auto bdofs_left = create_left_boundary_locator(V);

        for (auto &dof : bdofs_left) {
            dof = VxDofs.second[dof];
        }

        auto bdofs_right = create_right_boundary_locator(V);

        for (auto &dof : bdofs_right) {
            dof = VxDofs.second[dof];
        }

        auto bc1 = std::make_shared<const fem::DirichletBC<T>>(0, bdofs_left, X_sub);
        auto bc2 = std::make_shared<const fem::DirichletBC<T>>(imposed_disp[d], bdofs_right, X_sub);

        bcs.push_back(bc1);
        bcs.push_back(bc2);
    }

    // Phase BCs
    {
        auto CxDofs = X->sub({1})->collapse();
        auto &C = CxDofs.first;
        auto frac_locator = create_frac_locator(C);

        for (auto &dof : frac_locator) {
            dof = CxDofs.second[dof];
        }

        auto phase_field_BC = std::make_shared<const fem::DirichletBC<T>>(1., frac_locator, X->sub({1}));
        bcs.push_back(phase_field_BC);
    }

    return bcs;
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
        static int dims = 3;
        auto mesh =
            std::make_shared<mesh::Mesh<T>>(mesh::create_box(MPI_COMM_WORLD,
                                                             {{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}}},
                                                             {10, 10, 10},
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

        T disp_x = 0.001;
        int n_steps = 100;
        

        auto disp_bcs = disp_BC(V, disp_x/n_steps);
        auto phase_bcs = phase_BC(C);
        auto coupled_bcs = coupled_BC(X, disp_x/n_steps);

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

        // auto nls1 = std::make_shared<utopia::Newton<utopia::PetscMatrix>>(
        //     std::make_shared<utopia::Factorization<utopia::PetscMatrix, utopia::PetscVector>>());

        // auto nls2 = std::make_shared<utopia::Newton<utopia::PetscMatrix>>(
        //     std::make_shared<utopia::Factorization<utopia::PetscMatrix, utopia::PetscVector>>());

        // nls1->linear_solver()->verbose(true);
        // nls2->linear_solver()->verbose(true);

        // nls1->verbose(verbose);
        // nls2->verbose(verbose);

        auto f1_to_c12 = [](const utopia::PetscVector &in, utopia::PetscVector &out) {
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

        std::shared_ptr<utopia::NewtonInterface<utopia::PetscMatrix, utopia::PetscVector>> nlsolver;

        if  //
            // (true)  //
            (false)  //
        {
            auto ls = std::make_shared<utopia::BiCGStab<utopia::PetscMatrix, utopia::PetscVector, utopia::HOMEMADE>>();
            ls->max_it(5000);

            auto newton = std::make_shared<utopia::Newton<utopia::PetscMatrix>>(ls);

            auto params = utopia::param_list(   //
                utopia::param("damping", 0.5),  //
                utopia::param("max_it", 50),    //
                utopia::param("verbose", true)  //
            );

            newton->read(params);
            nlsolver = newton;

        } else if   //
            (true)  //
        // (false)  //
        {
            auto tfa = std::make_shared<utopia::TwoFieldAlternateMinimization<utopia::PetscMatrix>>(nls1, nls2);
            tfa->set_field_functions(f1, f2);
            tfa->set_transfers(f1_to_c12, f2_to_c12, c12_to_f1, c12_to_f2);
            tfa->verbose(verbose);
            tfa->max_it(400);
            // tfa->verbosity_level(utopia::VERBOSITY_LEVEL_DEBUG);
            // tfa->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);

            tfa->field1_diff_tol(1e-14);
            tfa->field2_diff_tol(1e-14);
            nlsolver = tfa;

        } else {
            auto lsc12 = std::make_shared<utopia::KSP_MF<utopia::PetscMatrix, utopia::PetscVector>>();
            lsc12->ksp_type("fgmres");
            lsc12->pc_type("none");

            lsc12->verbose(verbose);
            nls1->verbose(verbose);
            nls2->verbose(verbose);

            auto tfs = std::make_shared<utopia::TwoFieldSPIN<utopia::PetscMatrix>>(lsc12, nls1, nls2);
            tfs->set_field_functions(f1, f2);
            tfs->set_transfers(f1_to_c12, f2_to_c12, c12_to_f1, c12_to_f2);
            tfs->additive_precond(false);
            // tfs.additive_precond(true);
            // tfs.verbosity_level(utopia::VERBOSITY_LEVEL_DEBUG);
            nlsolver = tfs;
        }

        io::VTKFile file_u(mesh->comm(), "sub_u.pvd", "w");
        io::VTKFile file_c(mesh->comm(), "sub_c.pvd", "w");

        utopia::PrototypeFunction<utopia::PetscMatrix, utopia::PetscVector> f(
            [&](const utopia::PetscVector &x, T &value) -> bool { return c12->value(x, value); },
            [&](const utopia::PetscVector &x, utopia::PetscVector &g) -> bool { return c12->gradient(x, g); },
            [&](const utopia::PetscVector &x, utopia::PetscMatrix &H) -> bool { return c12->hessian(x, H); },
            [&](const T t) -> bool { 
                file_u.write<T>({*f1->u()}, t);
                file_c.write<T>({*f2->u()}, t);

                auto disp_bcs = disp_BC(V, t*disp_x/n_steps);
                auto phase_bcs = phase_BC(C);
                // auto coupled_bcs = coupled_BC(X, t*disp_x/n_steps);
                f1->set_boundary_conditions(disp_bcs);
                f2->set_boundary_conditions(phase_bcs);

                return true; },
            [=](utopia::PetscVector &x) { c12->create_vector(x); });

        utopia::PseudoTimeStepper<utopia::PetscMatrix, utopia::PetscVector> time_stepper(nlsolver);
        time_stepper.set_end_time(n_steps);
        time_stepper.solve(f, solution);


        // {
        //     // main problem
        //     auto u = c12->u()->sub({0}).collapse();
        //     io::VTKFile file_u(mesh->comm(), "u.pvd", "w");
        //     file_u.write<T>({u}, 0.0);

        //     auto c = c12->u()->sub({1}).collapse();
        //     io::VTKFile file_c(mesh->comm(), "c.pvd", "w");
        //     file_c.write<T>({c}, 0.0);
        // }
    }

    return utopia::Utopia::Finalize();
}
