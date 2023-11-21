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
#include "hyperelasticity.h"

#include "utopia.hpp"

#include "utopia_TwoFieldAlternateMinimization.hpp"
#include "utopia_TwoFieldSPIN.hpp"

using namespace dolfinx;
using T = PetscScalar;

class HyperElasticProblem : public utopia::Function<utopia::PetscMatrix, utopia::PetscVector> {
public:
    class Impl {
    public:
        Impl(std::shared_ptr<fem::FunctionSpace> V,
             std::shared_ptr<fem::Function<T>> u,
             std::shared_ptr<fem::Form<T>> objective,
             std::shared_ptr<fem::Form<T>> gradient,
             std::shared_ptr<fem::Form<T>> hessian,
             std::vector<std::shared_ptr<const fem::DirichletBC<T>>> boundary_conditions)
            : V(V),
              u(u),
              objective(objective),
              gradient(gradient),
              hessian(hessian),
              boundary_conditions(boundary_conditions),
              rhs(gradient->function_spaces()[0]->dofmap()->index_map,
                  gradient->function_spaces()[0]->dofmap()->index_map_bs()),
              matrix(la::petsc::Matrix(fem::petsc::create_matrix(*hessian, "baij"), false)),
              solution_vector(la::petsc::Vector(la::petsc::create_vector_wrap(*u->x()), false)) {
            auto map = gradient->function_spaces()[0]->dofmap()->index_map;
            const int bs = gradient->function_spaces()[0]->dofmap()->index_map_bs();
            std::int32_t size_local = bs * map->size_local();

            std::vector<PetscInt> ghosts(map->ghosts().begin(), map->ghosts().end());
            std::int64_t size_global = bs * map->size_global();
            VecCreateGhostBlockWithArray(map->comm(),
                                         bs,
                                         size_local,
                                         size_global,
                                         ghosts.size(),
                                         ghosts.data(),
                                         rhs.array().data(),
                                         &rhs_petsc_);
        }

        ~Impl() {
            if (rhs_petsc_) VecDestroy(&rhs_petsc_);
        }

        std::shared_ptr<fem::FunctionSpace> V;
        std::shared_ptr<fem::Function<T>> u;
        std::shared_ptr<fem::Form<T>> objective, gradient, hessian;
        std::vector<std::shared_ptr<const fem::DirichletBC<T>>> boundary_conditions;

        la::Vector<T> rhs;
        la::petsc::Matrix matrix;
        la::petsc::Vector solution_vector;

        Vec rhs_petsc_;

        auto form() {
            return [](Vec x) {
                VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
                VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);
            };
        }

        /// Compute F at current point x
        auto F() {
            return [&](const Vec x, Vec) {
                // Assemble b and update ghosts
                xtl::span<T> b(rhs.mutable_array());
                std::fill(b.begin(), b.end(), 0.0);
                fem::assemble_vector<T>(b, *gradient);
                VecGhostUpdateBegin(rhs_petsc_, ADD_VALUES, SCATTER_REVERSE);
                VecGhostUpdateEnd(rhs_petsc_, ADD_VALUES, SCATTER_REVERSE);

                // Set bcs
                Vec x_local;
                VecGhostGetLocalForm(x, &x_local);
                PetscInt n = 0;
                VecGetSize(x_local, &n);
                const T *array = nullptr;
                VecGetArrayRead(x_local, &array);
                fem::set_bc<T>(b, boundary_conditions, xtl::span<const T>(array, n), -1.0);
                VecRestoreArrayRead(x, &array);
            };
        }

        /// Compute J = F' at current point x
        auto J() {
            return [&](const Vec, Mat A) {
                MatZeroEntries(A);
                fem::assemble_matrix(la::petsc::Matrix::set_block_fn(A, ADD_VALUES), *hessian, boundary_conditions);
                MatAssemblyBegin(A, MAT_FLUSH_ASSEMBLY);
                MatAssemblyEnd(A, MAT_FLUSH_ASSEMBLY);
                fem::set_diagonal(
                    la::petsc::Matrix::set_fn(A, INSERT_VALUES), *hessian->function_spaces()[0], boundary_conditions);
                MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
                MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
            };
        }
    };

    HyperElasticProblem(const std::shared_ptr<fem::FunctionSpace> &V,
                        std::vector<std::shared_ptr<const fem::DirichletBC<T>>> boundary_conditions) {
        auto u = std::make_shared<fem::Function<T>>(V);

        auto objective = std::make_shared<fem::Form<T>>(
            fem::create_form<T>(*form_hyperelasticity_Pi, {}, {{"u", u}}, {}, {}, V->mesh()));

        auto gradient =
            std::make_shared<fem::Form<T>>(fem::create_form<T>(*form_hyperelasticity_F_form, {V}, {{"u", u}}, {}, {}));

        auto hessian = std::make_shared<fem::Form<T>>(
            fem::create_form<T>(*form_hyperelasticity_J_form, {V, V}, {{"u", u}}, {}, {}));

        impl_ = utopia::make_unique<Impl>(V, u, objective, gradient, hessian, boundary_conditions);
    }

    void create_vector(utopia::PetscVector &x) const override {
        Vec v = impl_->solution_vector.vec();
        x.wrap(v);
    }

    /// Destructor
    virtual ~HyperElasticProblem() {}

    bool hessian(const utopia::PetscVector &x, utopia::PetscMatrix &H) const override {
        UTOPIA_TRACE_SCOPE("HyperElasticProblem::hessian");

        if (x.raw_type() != impl_->solution_vector.vec()) {
            x.convert_to(impl_->solution_vector.vec());
        }

        auto H_fun = impl_->J();

        if (H.empty()) {
            auto m = impl_->matrix.mat();
            H.wrap(m);
        }

        impl_->form()(impl_->solution_vector.vec());

        H_fun(x.raw_type(), H.raw_type());

        // utopia::disp(H);
        return true;
    }

    bool value(const utopia::PetscVector &x, PetscScalar &value) const override {
        UTOPIA_TRACE_SCOPE("HyperElasticProblem::value");

        if (x.raw_type() != impl_->solution_vector.vec()) {
            x.convert_to(impl_->solution_vector.vec());
        }

        impl_->form()(impl_->solution_vector.vec());
        value = fem::assemble_scalar<T>(*impl_->objective);
        value = x.comm().sum(value);
        return true;
    }

    bool gradient(const utopia::PetscVector &x, utopia::PetscVector &g) const override {
        UTOPIA_TRACE_SCOPE("HyperElasticProblem::gradient");

        if (x.raw_type() != impl_->solution_vector.vec()) {
            x.convert_to(impl_->solution_vector.vec());
        }

        auto grad_fun = impl_->F();

        impl_->form()(impl_->solution_vector.vec());

        grad_fun(impl_->solution_vector.vec(), nullptr);

        la::petsc::Vector temp(la::petsc::create_vector_wrap(impl_->rhs), false);

        auto v = temp.vec();
        assert(v);
        g.copy_from(v);

        if (g.has_nan_or_inf()) {
            utopia::Utopia::Abort("HyperElasticProblem::gradient: Detected NaN!");
        }
        return true;
    }

    std::shared_ptr<fem::Function<T>> u() { return impl_->u; }

private:
    std::unique_ptr<Impl> impl_;
};

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
        auto mesh = std::make_shared<mesh::Mesh>(mesh::create_box(MPI_COMM_WORLD,
                                                                  {{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}}},
                                                                  {10, 10, 10},
                                                                  mesh::CellType::tetrahedron,
                                                                  mesh::GhostMode::none));

        auto V = std::make_shared<fem::FunctionSpace>(
            fem::create_functionspace(functionspace_form_hyperelasticity_F_form, "u", mesh));

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

        HyperElasticProblem problem(V, bcs);

        /////////////////////////////////////////////////////////
        /// Plug-in Utopia HERE

        bool verbose = true;

        // Split-fields
        auto f1 = std::make_shared<HyperElasticProblem>(V, bcs);
        auto f2 = std::make_shared<HyperElasticProblem>(V, bcs);

        // Global optimization problem (coupled fields)
        auto c12 = std::make_shared<HyperElasticProblem>(V, bcs);

        // auto nls1 = std::make_shared<utopia::Newton<utopia::PetscMatrix>>(
        //     std::make_shared<utopia::Factorization<utopia::PetscMatrix, utopia::PetscVector>>());

        // auto nls2 = std::make_shared<utopia::Newton<utopia::PetscMatrix>>(
        //     std::make_shared<utopia::Factorization<utopia::PetscMatrix, utopia::PetscVector>>());

        auto nls1 = std::make_shared<utopia::Newton<utopia::PetscMatrix>>(
            std::make_shared<utopia::ConjugateGradient<utopia::PetscMatrix, utopia::PetscVector, utopia::HOMEMADE>>());

        auto nls2 = std::make_shared<utopia::Newton<utopia::PetscMatrix>>(
            std::make_shared<utopia::ConjugateGradient<utopia::PetscMatrix, utopia::PetscVector, utopia::HOMEMADE>>());

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
            // auto lsc12 =
            // std::make_shared<utopia::ConjugateGradient<utopia::PetscMatrix, utopia::PetscVector,
            // utopia::HOMEMADE>>();
            // lsc12->apply_gradient_descent_step(true);
            // lsc12->set_preconditioner(std::make_shared<utopia::InvDiagPreconditioner<utopia::PetscMatrix,
            // utopia::PetscVector>>());

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

        // Compute Cauchy stress
        // Construct appropriate Basix element for stress
        constexpr auto family = basix::element::family::P;
        const auto cell_type = mesh::cell_type_to_basix_type(mesh->topology().cell_type());
        constexpr auto variant = basix::element::lagrange_variant::equispaced;
        constexpr int k = 0;
        constexpr bool discontinuous = true;

        const basix::FiniteElement S_element = basix::create_element(family, cell_type, k, variant, discontinuous);
        auto S = std::make_shared<fem::FunctionSpace>(
            fem::create_functionspace(mesh, S_element, pow(mesh->geometry().dim(), 2)));

        const auto sigma_expression =
            fem::create_expression<T>(*expression_hyperelasticity_sigma, {{"u", c12->u()}}, {}, mesh);

        auto sigma = fem::Function<T>(S);
        sigma.name = "cauchy_stress";
        sigma.interpolate(sigma_expression);

        // Save solution in VTK format
        io::VTKFile file_u(mesh->comm(), "u.pvd", "w");
        file_u.write<T>({*c12->u()}, 0.0);

        // Save Cauchy stress in XDMF format
        io::XDMFFile file_sigma(mesh->comm(), "sigma.xdmf", "w");
        file_sigma.write_mesh(*mesh);
        file_sigma.write_function(sigma, 0.0);
    }

    // PetscFinalize();
    return utopia::Utopia::Finalize();
}
