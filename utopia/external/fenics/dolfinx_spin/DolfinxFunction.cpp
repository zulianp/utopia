#include "DolfinxFunction.hpp"

using namespace dolfinx;
using T = PetscScalar;

class DolfinxFunction::Impl {
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
        VecCreateGhostBlockWithArray(
            map->comm(), bs, size_local, size_global, ghosts.size(), ghosts.data(), rhs.array().data(), &rhs_petsc_);
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

    void set_BC(Vec x)
    {
        // Set bcs
        // Vec x_local;
        // VecGhostGetLocalForm(x, &x_local);
        // PetscInt n = 0;
        // VecGetSize(x_local, &n);
        // T *array = nullptr;
        // VecGetArray(x_local, &array);

        // xtl::span<T> b(array, n);

        // // FIXME Boundary conditions should be set to if x satisfies them (!!!)
        // fem::set_bc<T>(b, boundary_conditions, xtl::span<const T>(array, n), 1.0);
        // VecRestoreArray(x, &array);
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

            // FIXME Boundary conditions should be set to if x satisfies them (!!!)
            fem::set_bc<T>(b, boundary_conditions, xtl::span<const T>(array, n), -1.0);
            // fem::set_bc<T>(b, boundary_conditions, xtl::span<const T>(array, n), 1.0);
            // fem::set_bc<T>(b, boundary_conditions, xtl::span<const T>(array, n), 0);
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

DolfinxFunction::DolfinxFunction(std::shared_ptr<fem::FunctionSpace> V,
                                 std::shared_ptr<fem::Function<T>> u,
                                 std::shared_ptr<fem::Form<T>> objective,
                                 std::shared_ptr<fem::Form<T>> gradient,
                                 std::shared_ptr<fem::Form<T>> hessian,
                                 std::vector<std::shared_ptr<const fem::DirichletBC<T>>> boundary_conditions) {
    impl_ = utopia::make_unique<Impl>(V, u, objective, gradient, hessian, boundary_conditions);
}

DolfinxFunction::~DolfinxFunction() {}

void DolfinxFunction::create_vector(utopia::PetscVector &x) const {
    Vec v = impl_->solution_vector.vec();

    impl_->set_BC(v);
    x.wrap(v);

    utopia::out() << " DolfinxFunction::create_vector " << x.size() << "\n";
}

bool DolfinxFunction::hessian(const utopia::PetscVector &x, utopia::PetscMatrix &H) const {
    UTOPIA_TRACE_SCOPE("DolfinxFunction::hessian");

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
    return true;
}

bool DolfinxFunction::value(const utopia::PetscVector &x, PetscScalar &value) const {
    UTOPIA_TRACE_SCOPE("DolfinxFunction::value");

    if (x.raw_type() != impl_->solution_vector.vec()) {
        x.convert_to(impl_->solution_vector.vec());
    }

    impl_->form()(impl_->solution_vector.vec());
    value = fem::assemble_scalar<T>(*impl_->objective);
    value = x.comm().sum(value);
    return true;
}

bool DolfinxFunction::gradient(const utopia::PetscVector &x, utopia::PetscVector &g) const {
    UTOPIA_TRACE_SCOPE("DolfinxFunction::gradient");

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
        utopia::Utopia::Abort("DolfinxFunction::gradient: Detected NaN!");
    }
    return true;
}

std::shared_ptr<fem::Function<T>> DolfinxFunction::u() { return impl_->u; }
