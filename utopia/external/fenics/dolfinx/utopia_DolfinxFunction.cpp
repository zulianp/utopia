#include "utopia_DolfinxFunction.hpp"

#include "utopia_make_unique.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Vector.hpp"

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

using namespace dolfinx;

namespace utopia {
    class DolfinxFunction::Impl {
    public:
        std::shared_ptr<Form_t> value;
        std::shared_ptr<Form_t> gradient;
        std::shared_ptr<Form_t> hessian;
        std::vector<std::shared_ptr<const DirichletBC_t>> bcs;
        la::Vector<Scalar> gradient_vector;
        Vec gradient_vector_ghosts = nullptr;
        la::petsc::Matrix hessian_matrix;

        Impl(
        	std::shared_ptr<Form_t> value,
        	std::shared_ptr<Form_t> gradient,
        	std::shared_ptr<Form_t> hessian,
        	std::vector<std::shared_ptr<const DirichletBC_t>> bcs)
        : value(value), gradient(gradient), hessian(hessian), bcs(bcs),  gradient_vector(gradient->function_spaces()[0]->dofmap()->index_map, gradient->function_spaces()[0]->dofmap()->index_map_bs()),
          hessian_matrix(la::petsc::Matrix(fem::petsc::create_matrix(*hessian, "baij"), false))
          {
          	auto map = gradient->function_spaces()[0]->dofmap()->index_map;
          	const int bs = gradient->function_spaces()[0]->dofmap()->index_map_bs();
          	std::int32_t size_local = bs * map->size_local();

          	std::vector<PetscInt> ghosts(map->ghosts().begin(), map->ghosts().end());
          	std::int64_t size_global = bs * map->size_global();
          	VecCreateGhostBlockWithArray(
          	    map->comm(), bs, size_local, size_global, ghosts.size(), ghosts.data(), gradient_vector.array().data(), &gradient_vector_ghosts);
          }

        void update_ghosts(Vec x) {
            VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
            VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);
        }

        void eval_gradient(const Vec x) {
            xtl::span<Scalar> b(gradient_vector.mutable_array());
            std::fill(b.begin(), b.end(), 0.0);
            fem::assemble_vector<Scalar>(b, *gradient);
            VecGhostUpdateBegin(gradient_vector_ghosts, ADD_VALUES, SCATTER_REVERSE);
            VecGhostUpdateEnd(gradient_vector_ghosts, ADD_VALUES, SCATTER_REVERSE);

            // Set bcs
            Vec x_local;
            VecGhostGetLocalForm(x, &x_local);
            PetscInt n = 0;
            VecGetSize(x_local, &n);
            const Scalar *array = nullptr;
            VecGetArrayRead(x_local, &array);
            fem::set_bc<Scalar>(b, bcs, xtl::span<const Scalar>(array, n), -1.0);
            VecRestoreArrayRead(x, &array);
        }

        //FIXME! x is not used!
        void eval_hessian(const Vec x, Mat A) {
            MatZeroEntries(A);
            fem::assemble_matrix(la::petsc::Matrix::set_block_fn(A, ADD_VALUES), *hessian, bcs);
            MatAssemblyBegin(A, MAT_FLUSH_ASSEMBLY);
            MatAssemblyEnd(A, MAT_FLUSH_ASSEMBLY);
            fem::set_diagonal(la::petsc::Matrix::set_fn(A, INSERT_VALUES), *hessian->function_spaces()[0], bcs);
            MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
        }

        bool eval_value(const Vec x, PetscScalar &value) const  {
            value = fem::assemble_scalar<Scalar>(*this->value);
            return true;
        }
    };

    DolfinxFunction::DolfinxFunction(const std::shared_ptr<Form_t> &value,
                                     const std::shared_ptr<Form_t> &gradient,
                                     const std::shared_ptr<Form_t> &hessian,
                                     std::vector<std::shared_ptr<const DirichletBC_t>> bcs)
        : impl_(utopia::make_unique<Impl>(value, gradient, hessian, bcs)) {
    }

    DolfinxFunction::DolfinxFunction(const std::shared_ptr<Form_t> &residual,
                                     const std::shared_ptr<Form_t> &jacobian,
                                     std::vector<std::shared_ptr<const DirichletBC_t>> bcs)
        : impl_(utopia::make_unique<Impl>(nullptr, residual, jacobian, bcs)) {
    }

    DolfinxFunction::~DolfinxFunction() = default;

    void DolfinxFunction::init() {
        // TODO
    }

    void DolfinxFunction::destory() {
        if (impl_->gradient_vector_ghosts) VecDestroy(&impl_->gradient_vector_ghosts);
    }

    bool DolfinxFunction::hessian(const utopia::PetscVector &x, utopia::PetscMatrix &H) const {
        if (H.empty()) {
        	auto m = impl_->hessian_matrix.mat();
            H.wrap(m);
        }

        impl_->update_ghosts(x.raw_type());
        impl_->eval_hessian(x.raw_type(), H.raw_type());
        return true;
    }

    bool DolfinxFunction::value(const utopia::PetscVector &x, PetscScalar &value) const {
        impl_->update_ghosts(x.raw_type());
        impl_->eval_value(x.raw_type(), value);
        value = x.comm().sum(value);
        return true;
    }

    bool DolfinxFunction::gradient(const utopia::PetscVector &x, utopia::PetscVector &g) const {
        impl_->update_ghosts(x.raw_type());
        impl_->eval_gradient(x.raw_type());

        la::petsc::Vector temp(la::petsc::create_vector_wrap(impl_->gradient_vector), false);

        auto v = temp.vec();
        assert(v);
        g.copy_from(v);
        return true;
    }

}  // namespace utopia
