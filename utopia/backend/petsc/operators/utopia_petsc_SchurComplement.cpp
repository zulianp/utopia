#include "utopia_petsc_SchurComplement.hpp"

#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Vector.hpp"

#include "utopia_DeviceView.hpp"

#include "utopia_Instance.hpp"
#include "utopia_petsc_Factorization.hpp"

namespace utopia {

    class SchurComplement<PetscMatrix>::Impl {
    public:
        IS is_eliminated{nullptr};
        IS is_dof{nullptr};

        PetscMatrix A_II, A_IG, A_GI, A_GG;
        PetscVector x_I, x_G, temp_I;
        Mat *submat{nullptr};

        std::shared_ptr<LinearSolver<PetscMatrix, PetscVector>> A_II_inv;

        void destroy() {
            if (is_eliminated) ISDestroy(&is_eliminated);
            if (is_dof) ISDestroy(&is_dof);
            A_II_inv = nullptr;

            destroy_matrices();
        }

        void destroy_matrices() {
            if (submat) {
                MatDestroySubMatrices(4, &submat);
            }
        }

        bool create_A_II_inverse() {
            UTOPIA_TRACE_REGION_BEGIN("SchurComplement::initialize::decomposition");
            if (!A_II_inv) {
                A_II_inv = std::make_shared<Factorization<PetscMatrix, PetscVector>>("mumps", "cholesky");
            }

            A_II_inv->update(make_ref(A_II));
            UTOPIA_TRACE_REGION_END("SchurComplement::initialize::decomposition");
            return true;
        }

        bool create_index_sets(const PetscMatrix &matrix, const IndexArray &eliminated_dofs) {
            this->destroy();

            PetscErrorCode err = ISCreateGeneral(matrix.comm().raw_comm(),
                                                 eliminated_dofs.size(),
                                                 &eliminated_dofs[0],
                                                 PETSC_USE_POINTER,
                                                 &this->is_eliminated);
            if (err != 0) {
                return false;
            }

            std::vector<bool> selected(matrix.local_rows(), false);
            auto rr = matrix.row_range();

            for (auto e : eliminated_dofs) {
                selected[e - rr.begin()] = true;
            }

            size_t n_free = rr.extent() - eliminated_dofs.size();

            PetscInt *free_dofs = nullptr;
            PetscMalloc1(n_free * sizeof(PetscInt), &free_dofs);

            PetscInt n = matrix.local_rows();
            PetscInt offset = 0;

            for (PetscInt i = 0; i < n; ++i) {
                if (selected[i]) {
                    free_dofs[offset++] = i + rr.begin();
                }
            }

            err = ISCreateGeneral(
                matrix.comm().raw_comm(), eliminated_dofs.size(), free_dofs, PETSC_OWN_POINTER, &this->is_dof);

            return true;
        }

        bool create_matrix_blocks(const PetscMatrix &matrix, MatReuse scall = MAT_INITIAL_MATRIX) {
            IS irow[4] = {is_eliminated, is_eliminated, is_dof, is_dof};
            IS icol[4] = {is_eliminated, is_dof, is_eliminated, is_dof};

            PetscErrorCode err;

            if (scall == MAT_INITIAL_MATRIX) {
                err = MatCreateSubMatrices(matrix.raw_type(), 4, irow, icol, scall, &submat);
            } else {
                assert(false && "IMPLEMENT ME");
                Utopia::Abort("Not implemented!");
            }

            A_II.wrap(submat[0]);
            A_IG.wrap(submat[1]);
            A_GI.wrap(submat[2]);
            A_GG.wrap(submat[3]);
            return true;
        }

        bool create_vector_block(const PetscVector &vector, IS is, PetscVector &subv) {
            assert(empty(subv) || subv.comm().size() == 1);

            PetscInt size_is;
            ISGetLocalSize(is, &size_is);

            if (empty(subv) || subv.local_size() != size_is) {
                subv.zeros(serial_layout(size_is));
            }

            auto v_view = const_local_view_device(vector);
            auto subv_view = local_view_device(subv);

            const PetscInt *idx = nullptr;
            ISGetIndices(is, &idx);

            parallel_for(
                local_range_device(subv), UTOPIA_LAMBDA(const SizeType i) { subv_view.set(i, v_view.get(idx[i])); });

            ISRestoreIndices(is, &idx);
            return true;
        }

        bool create_vector_blocks(const PetscVector &vector, PetscVector &v_I, PetscVector &v_G) {
            return create_vector_block(vector, is_eliminated, v_I) && create_vector_block(vector, is_dof, v_G);
        }

        bool restore_vector_blocks(PetscVector &vector, const PetscVector &v_I, const PetscVector &v_G) {
            assert(!empty(vector));
            assert(vector.local_size() == v_I.local_size() + v_G.local_size());

            auto v_view = local_view_device(vector);
            auto v_I_view = const_local_view_device(v_I);
            auto v_G_view = const_local_view_device(v_G);

            const PetscInt *idx_dof = nullptr;
            ISGetIndices(is_dof, &idx_dof);

            const PetscInt *idx_eliminated = nullptr;
            ISGetIndices(is_eliminated, &idx_eliminated);

            parallel_for(
                local_range_device(v_I),
                UTOPIA_LAMBDA(const SizeType i) { v_view.set(idx_eliminated[i], v_I_view.get(i)); });

            parallel_for(
                local_range_device(v_G), UTOPIA_LAMBDA(const SizeType i) { v_view.set(idx_dof[i], v_G_view.get(i)); });

            ISRestoreIndices(is_eliminated, &idx_eliminated);
            ISRestoreIndices(is_dof, &idx_dof);
            return true;
        }
    };

    void SchurComplement<PetscMatrix>::read(Input &) {
        // TODO
    }

    SchurComplement<PetscMatrix>::SchurComplement() : impl_(utopia::make_unique<Impl>()) {}
    SchurComplement<PetscMatrix>::~SchurComplement() = default;

    bool SchurComplement<PetscMatrix>::apply_righthand_side(const PetscVector &rhs, PetscVector &out) {
        if (empty(impl_->temp_I) || impl_->temp_I.local_size() != impl_->A_GG.local_rows()) {
            impl_->temp_I.zeros(row_layout(impl_->A_GG));
        } else {
            impl_->temp_I.set(0.0);
        }

        impl_->create_vector_blocks(rhs, impl_->x_G, impl_->x_I);

        if (!impl_->A_II_inv->apply(impl_->x_I, impl_->temp_I)) {
            // return false;
        }

        out = impl_->A_GI * impl_->temp_I;
        out = impl_->x_G - out;
        return true;
    }

    bool SchurComplement<PetscMatrix>::initialize_from_selection(const PetscMatrix &matrix,
                                                                 const IndexArray &eliminated_dofs) {
        if (!impl_->create_index_sets(matrix, eliminated_dofs)) {
            return false;
        }

        return impl_->create_matrix_blocks(matrix) && impl_->create_A_II_inverse();
    }

    void SchurComplement<PetscMatrix>::set_solver(
        const std::shared_ptr<LinearSolver<PetscMatrix, PetscVector>> &solver) {
        impl_->A_II_inv = solver;
    }

    bool SchurComplement<PetscMatrix>::apply(const PetscVector &x, PetscVector &y) const {
        impl_->create_vector_block(x, impl_->is_dof, impl_->x_G);
        impl_->temp_I = impl_->A_IG * impl_->x_G;

        if (empty(impl_->x_I) || impl_->x_I.local_size() != impl_->temp_I.local_size()) {
            impl_->x_I.zeros(layout(impl_->temp_I));
        } else {
            impl_->x_I.set(0.0);
        }

        if (impl_->A_II_inv->apply(impl_->temp_I, impl_->x_I)) {
            return false;
        }

        y = impl_->A_GG * impl_->x_G;
        y -= impl_->A_GI * impl_->x_I;
        return true;
    }

}  // namespace utopia
