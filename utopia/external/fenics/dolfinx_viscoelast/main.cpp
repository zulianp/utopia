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

#include "DolfinxFunction.hpp"

#include <climits>
#include <cmath>
#include "NLViscoElast.h"

#include "utopia.hpp"

using namespace dolfinx;
using T = PetscScalar;

std::shared_ptr<DolfinxFunction> NLViscoElast(
    const std::shared_ptr<fem::FunctionSpace<T>> &V,
    std::shared_ptr<dolfinx::fem::Function<T>> u,
    std::shared_ptr<dolfinx::fem::Function<T>> u_old,
    std::vector<std::shared_ptr<const fem::DirichletBC<T>>> boundary_conditions) {
    auto objective = std::make_shared<fem::Form<T>>(
        fem::create_form<T>(*form_NLViscoElast_energy, {}, {{"u", u}, {"u_old", u_old}}, {}, {}, V->mesh()));

    auto gradient = std::make_shared<fem::Form<T>>(
        fem::create_form<T>(*form_NLViscoElast_gradient, {V}, {{"u", u}, {"u_old", u_old}}, {}, {}));

    auto hessian = std::make_shared<fem::Form<T>>(
        fem::create_form<T>(*form_NLViscoElast_hessian, {V, V}, {{"u", u}, {"u_old", u_old}}, {}, {}));

    return std::make_shared<DolfinxFunction>(V, u, objective, gradient, hessian, boundary_conditions);
}

int main(int argc, char *argv[]) {
    dolfinx::init_logging(argc, argv);
    // PetscInitialize(&argc, &argv, nullptr, nullptr);
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
        auto mesh =
            std::make_shared<mesh::Mesh<T>>(mesh::create_box(MPI_COMM_WORLD,
                                                             {{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}}},
                                                             {10, 10, 10},
                                                             mesh::CellType::tetrahedron,
                                                             mesh::create_cell_partitioner(mesh::GhostMode::none)));

        auto V = std::make_shared<fem::FunctionSpace<T>>(
            fem::create_functionspace(functionspace_form_NLViscoElast_gradient, "u", mesh));

        // Define solution function
        auto u = std::make_shared<fem::Function<T>>(V);
        auto u_old = std::make_shared<fem::Function<T>>(V);


        // auto u_rotation = std::make_shared<fem::Function<T>>(V);
        // u_rotation->interpolate([](auto x) -> std::pair<std::vector<T>, std::vector<std::size_t>> {
        //     constexpr T scale = 0.005;

        //     // Center of rotation
        //     constexpr T x1_c = 0.5;
        //     constexpr T x2_c = 0.5;

        //     // Large angle of rotation (60 degrees)
        //     constexpr T theta = 1.04719755;

        //     // New coordinates
        //     std::vector<T> fdata(3 * x.extent(1), 0.0);
        //     MDSPAN_IMPL_STANDARD_NAMESPACE::mdspan<
        //         T,
        //         MDSPAN_IMPL_STANDARD_NAMESPACE::extents<std::size_t, 3, MDSPAN_IMPL_STANDARD_NAMESPACE::dynamic_extent>>
        //         f(fdata.data(), 3, x.extent(1));
        //     for (std::size_t p = 0; p < x.extent(1); ++p) {
        //         T x1 = x(1, p);
        //         T x2 = x(2, p);
        //         f(1, p) = scale * (x1_c + (x1 - x1_c) * std::cos(theta) - (x2 - x2_c) * std::sin(theta) - x1);
        //         f(2, p) = scale * (x2_c + (x1 - x1_c) * std::sin(theta) - (x2 - x2_c) * std::cos(theta) - x2);
        //     }

        //     return {std::move(fdata), {3, x.extent(1)}};
        // });

        // Create Dirichlet boundary conditions
        auto bdofs_left = fem::locate_dofs_geometrical(*V, [](auto x) {
            constexpr T eps = 1.0e-6;
            std::vector<std::int8_t> marker(x.extent(1), false);
            for (std::size_t p = 0; p < x.extent(1); ++p) {
                if (std::abs(x(0, p)) < eps) marker[p] = true;
            }
            return marker;
        });
        auto bdofs_right = fem::locate_dofs_geometrical(*V, [](auto x) {
            constexpr T eps = 1.0e-6;
            std::vector<std::int8_t> marker(x.extent(1), false);
            for (std::size_t p = 0; p < x.extent(1); ++p) {
                if (std::abs(x(0, p) - 1) < eps) marker[p] = true;
            }
            return marker;
        });
        auto bcs = std::vector{std::make_shared<const fem::DirichletBC<T>>(std::vector<T>{0, 0, 0}, bdofs_left, V),
                               // std::make_shared<const fem::DirichletBC<T>>(u_rotation, bdofs_right)
                               std::make_shared<const fem::DirichletBC<T>>(std::vector<T>{0.5, 0, 0}, bdofs_right, V)
                           };

        auto problem = NLViscoElast(V, u, u_old, bcs);

        utopia::Newton<utopia::PetscMatrix, utopia::PetscVector> newton(
            std::make_shared<utopia::ConjugateGradient<utopia::PetscMatrix, utopia::PetscVector, utopia::HOMEMADE>>()
            // std::make_shared<utopia::BiCGStab<utopia::PetscMatrix, utopia::PetscVector>>()
        );

        newton.verbose(true);
        utopia::PetscVector solution;
        problem->create_vector(solution);

        // io::VTKFile file_u(mesh->comm(), "steady.pvd", "w");
        io::VTKFile file_u(mesh->comm(), "u.pvd", "w");

        Vec u_petsc = la::petsc::create_vector_wrap(*u->x());
        Vec u_old_petsc = la::petsc::create_vector_wrap(*u_old->x());

        for(int t = 0; t < 100; t++) {
            newton.solve(*problem, solution);
            file_u.write<T>({*u}, t * 0.01);

            VecCopy(u_petsc, u_old_petsc);

        }
    }

    // PetscFinalize();
    return utopia::Utopia::Finalize();
}
