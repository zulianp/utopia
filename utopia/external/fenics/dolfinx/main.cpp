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

using namespace dolfinx;
using T = PetscScalar;

class HyperElasticProblem : public utopia::Function<utopia::PetscMatrix, utopia::PetscVector> {
public:
    HyperElasticProblem(std::shared_ptr<fem::Form<T>> e,
                        std::shared_ptr<fem::Form<T>> L,
                        std::shared_ptr<fem::Form<T>> J,
                        std::vector<std::shared_ptr<const fem::DirichletBC<T>>> bcs)
        : _e(e),
          _l(L),
          _j(J),
          _bcs(bcs),
          _b(L->function_spaces()[0]->dofmap()->index_map, L->function_spaces()[0]->dofmap()->index_map_bs()),
          _matA(la::petsc::Matrix(fem::petsc::create_matrix(*J, "baij"), false)) {
        auto map = L->function_spaces()[0]->dofmap()->index_map;
        const int bs = L->function_spaces()[0]->dofmap()->index_map_bs();
        std::int32_t size_local = bs * map->size_local();

        std::vector<PetscInt> ghosts(map->ghosts().begin(), map->ghosts().end());
        std::int64_t size_global = bs * map->size_global();
        VecCreateGhostBlockWithArray(
            map->comm(), bs, size_local, size_global, ghosts.size(), ghosts.data(), _b.array().data(), &_b_petsc);
    }

    /// Destructor
    virtual ~HyperElasticProblem() {
        if (_b_petsc) VecDestroy(&_b_petsc);
    }

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
            xtl::span<T> b(_b.mutable_array());
            std::fill(b.begin(), b.end(), 0.0);
            fem::assemble_vector<T>(b, *_l);
            VecGhostUpdateBegin(_b_petsc, ADD_VALUES, SCATTER_REVERSE);
            VecGhostUpdateEnd(_b_petsc, ADD_VALUES, SCATTER_REVERSE);

            // Set bcs
            Vec x_local;
            VecGhostGetLocalForm(x, &x_local);
            PetscInt n = 0;
            VecGetSize(x_local, &n);
            const T *array = nullptr;
            VecGetArrayRead(x_local, &array);
            fem::set_bc<T>(b, _bcs, xtl::span<const T>(array, n), -1.0);
            VecRestoreArrayRead(x, &array);
        };
    }

    /// Compute J = F' at current point x
    auto J() {
        return [&](const Vec, Mat A) {
            MatZeroEntries(A);
            fem::assemble_matrix(la::petsc::Matrix::set_block_fn(A, ADD_VALUES), *_j, _bcs);
            MatAssemblyBegin(A, MAT_FLUSH_ASSEMBLY);
            MatAssemblyEnd(A, MAT_FLUSH_ASSEMBLY);
            fem::set_diagonal(la::petsc::Matrix::set_fn(A, INSERT_VALUES), *_j->function_spaces()[0], _bcs);
            MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
        };
    }

    bool hessian(const utopia::PetscVector &x, utopia::PetscMatrix &H) const override {
        auto that = const_cast<HyperElasticProblem *>(this);
        auto H_fun = that->J();

        if (H.empty()) {
            auto m = that->matrix();
            H.wrap(m);
        }

        that->form()(x.raw_type());

        H_fun(x.raw_type(), H.raw_type());
        return true;
    }

    bool value(const utopia::PetscVector &x, PetscScalar &value) const override {
        auto that = const_cast<HyperElasticProblem *>(this);
        that->form()(x.raw_type());
        value = fem::assemble_scalar<T>(*_e);
        value = x.comm().sum(value);
        return true;
    }

    bool gradient(const utopia::PetscVector &x, utopia::PetscVector &g) const override {
        auto that = const_cast<HyperElasticProblem *>(this);
        auto grad_fun = that->F();

        that->form()(x.raw_type());

        grad_fun(x.raw_type(), nullptr);

        la::petsc::Vector temp(la::petsc::create_vector_wrap(_b), false);

        auto v = temp.vec();
        assert(v);
        g.copy_from(v);
        return true;
    }

    Vec vector() const { return _b_petsc; }

    Mat matrix() const { return _matA.mat(); }

private:
    std::shared_ptr<fem::Form<T>> _e, _l, _j;
    std::vector<std::shared_ptr<const fem::DirichletBC<T>>> _bcs;
    la::Vector<T> _b;
    Vec _b_petsc = nullptr;
    la::petsc::Matrix _matA;
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

        // Define solution function
        auto u = std::make_shared<fem::Function<T>>(V);

        auto e = std::make_shared<fem::Form<T>>(
            fem::create_form<T>(*form_hyperelasticity_Pi, {}, {{"u", u}}, {}, {}, mesh)
            );
        
        auto a = std::make_shared<fem::Form<T>>(
            fem::create_form<T>(*form_hyperelasticity_J_form, {V, V}, {{"u", u}}, {}, {}));
        
        auto L =
            std::make_shared<fem::Form<T>>(fem::create_form<T>(*form_hyperelasticity_F_form, {V}, {{"u", u}}, {}, {}));

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

        HyperElasticProblem problem(e, L, a, bcs);

        la::petsc::Vector _u(la::petsc::create_vector_wrap(*u->x()), false);

        /////////////////////////////////////////////////////////
        /// Plug-in Utopia HERE

        if (false) {
            nls::petsc::NewtonSolver newton_solver(mesh->comm());
            newton_solver.setF(problem.F(), problem.vector());
            newton_solver.setJ(problem.J(), problem.matrix());
            newton_solver.set_form(problem.form());
            newton_solver.solve(_u.vec());
        } else {
            utopia::Newton<utopia::PetscMatrix, utopia::PetscVector> newton;
            newton.verbose(true);
            // auto params = utopia::param_list(utopia::param("damping", 0.5));
            // newton.read(params);

            utopia::PetscVector uvec;
            {
                auto u_ = _u.vec();
                uvec.wrap(u_);
            }

            newton.solve(problem, uvec);
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
            fem::create_expression<T>(*expression_hyperelasticity_sigma, {{"u", u}}, {}, mesh);

        auto sigma = fem::Function<T>(S);
        sigma.name = "cauchy_stress";
        sigma.interpolate(sigma_expression);

        // Save solution in VTK format
        io::VTKFile file_u(mesh->comm(), "u.pvd", "w");
        file_u.write<T>({*u}, 0.0);

        // Save Cauchy stress in XDMF format
        io::XDMFFile file_sigma(mesh->comm(), "sigma.xdmf", "w");
        file_sigma.write_mesh(*mesh);
        file_sigma.write_function(sigma, 0.0);
    }

    // PetscFinalize();
    return utopia::Utopia::Finalize();
}
