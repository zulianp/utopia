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
        PetscVector x_I_local, temp_I_local;

        Mat *submat{nullptr};

        std::shared_ptr<LinearSolver<PetscMatrix, PetscVector>> A_II_inv;

        Communicator comm;

        static void copy_values(const PetscVector &x_from, PetscVector &x_to) {
            {
                // FIXME (copying stuff just because of abstractions!)
                auto x_from_view = local_view_device(x_from);
                auto x_to_view = local_view_device(x_to);
                parallel_for(
                    local_range_device(x_to),
                    UTOPIA_LAMBDA(const SizeType i) { x_to_view.set(i, x_from_view.get(i)); });
            }
        }

        static void parallel_to_serial(const PetscVector &x_from, PetscVector &x_to) {
            if (empty(x_to)) {
                x_to.zeros(serial_layout(x_from.local_size()));
            }

            copy_values(x_from, x_to);
        }

        static void serial_to_parallel(const PetscVector &x_from, PetscVector &x_to) { copy_values(x_from, x_to); }

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

            assert(A_II.comm().size() == 1);
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

            if (n_free) {
                for (PetscInt i = 0; i < n; ++i) {
                    if (!selected[i]) {
                        free_dofs[offset++] = i + rr.begin();
                    }
                }
            }

            assert(size_t(offset) == n_free);

            err = ISCreateGeneral(matrix.comm().raw_comm(), n_free, free_dofs, PETSC_OWN_POINTER, &this->is_dof);

            return true;
        }

        bool create_matrix_blocks(const PetscMatrix &matrix, MatReuse scall = MAT_INITIAL_MATRIX) {
            IS irow[4] = {is_eliminated, is_eliminated, is_dof, is_dof};
            IS icol[4] = {is_eliminated, is_dof, is_eliminated, is_dof};

            PetscErrorCode err = 0;

            this->comm = matrix.comm();

            if (matrix.comm().size() == 1) {
                // if (false) {
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
            } else {
                if (scall == MAT_INITIAL_MATRIX) {
                    A_IG.destroy();
                    MatCreateSubMatrix(matrix.raw_type(), is_eliminated, is_dof, scall, &A_IG.raw_type());

                    A_GI.destroy();
                    MatCreateSubMatrix(matrix.raw_type(), is_dof, is_eliminated, scall, &A_GI.raw_type());

                    A_GG.destroy();
                    MatCreateSubMatrix(matrix.raw_type(), is_dof, is_dof, scall, &A_GG.raw_type());

                    // Serial matrix!
                    Mat *temp{nullptr};
                    MatCreateSubMatrices(matrix.raw_type(), 1, &is_eliminated, &is_eliminated, scall, &temp);

                    A_II.own(temp[0]);
                    PetscFree(temp);
                } else {
                    assert(false && "IMPLEMENT ME");
                    Utopia::Abort("Not implemented!");
                }
            }

            return err == 0;
        }

        bool create_parallel_vector_block(const PetscVector &vector, IS is, PetscVector &subv) {
            PetscInt size_is;
            ISGetLocalSize(is, &size_is);

            if (empty(subv) || vector.comm().disjunction(subv.local_size() != size_is)) {
                subv.zeros(layout(vector.comm(), size_is, Traits<PetscMatrix>::determine()));
            }

            auto v_view = const_local_view_device(vector);
            auto subv_view = local_view_device(subv);

            const PetscInt *idx = nullptr;
            ISGetIndices(is, &idx);

            PetscInt offset = range(vector).begin();
            parallel_for(
                local_range_device(subv),
                UTOPIA_LAMBDA(const SizeType i) { subv_view.set(i, v_view.get(idx[i] - offset)); });

            ISRestoreIndices(is, &idx);
            return true;
        }

        bool create_vector_blocks(const PetscVector &vector, PetscVector &v_I, PetscVector &v_G) {
            return create_parallel_vector_block(vector, is_eliminated, v_I) &&
                   create_parallel_vector_block(vector, is_dof, v_G);
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

            PetscInt offset = range(vector).begin();

            parallel_for(
                local_range_device(v_I),
                UTOPIA_LAMBDA(const SizeType i) { v_view.set(idx_eliminated[i] - offset, v_I_view.get(i)); });

            parallel_for(
                local_range_device(v_G),
                UTOPIA_LAMBDA(const SizeType i) { v_view.set(idx_dof[i] - offset, v_G_view.get(i)); });

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

    bool SchurComplement<PetscMatrix>::apply_righthand_side(const PetscVector &rhs,
                                                            PetscVector &out_eliminated,
                                                            PetscVector &out_restricted) {
        UTOPIA_TRACE_REGION_BEGIN("SchurComplement::apply_righthand_side");

        if (empty(impl_->temp_I) || !layout(impl_->temp_I).same(row_layout(impl_->A_II))) {
            impl_->temp_I.zeros(row_layout(impl_->A_IG));
        } else {
            impl_->temp_I.set(0.0);
        }

        impl_->create_vector_blocks(rhs, out_eliminated, impl_->x_G);

        impl_->temp_I.create_local_vector(impl_->temp_I_local);
        out_eliminated.create_local_vector(impl_->x_I_local);

        assert(impl_->x_I_local.comm().size() == 1);
        assert(impl_->temp_I_local.comm().size() == 1);

        if (!impl_->A_II_inv->apply(impl_->x_I_local, impl_->temp_I_local)) {
            // return false;
        }

        impl_->temp_I.restore_local_vector(impl_->temp_I_local);
        out_eliminated.restore_local_vector(impl_->x_I_local);

        if (out_restricted.is_alias(rhs)) {
            rhs.comm().root_print("Alias!", utopia::out().stream());
        }

        if (empty(out_restricted) || !layout(out_restricted).same(row_layout(impl_->A_GI))) {
            out_restricted.zeros(row_layout(impl_->A_GI));
        }

        out_restricted = impl_->A_GI * impl_->temp_I;
        out_restricted = impl_->x_G - out_restricted;

        UTOPIA_TRACE_REGION_END("SchurComplement::apply_righthand_side");
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
        impl_->x_G = x;
        impl_->temp_I = impl_->A_IG * impl_->x_G;

        if (empty(impl_->x_I) || !layout(impl_->x_I).same(layout(impl_->temp_I))) {
            impl_->x_I.zeros(layout(impl_->temp_I));
        } else {
            impl_->x_I.set(0.0);
        }

        impl_->temp_I.create_local_vector(impl_->temp_I_local);
        impl_->x_I.create_local_vector(impl_->x_I_local);

        assert(impl_->temp_I_local.comm().size() == 1);
        assert(impl_->x_I_local.comm().size() == 1);

        if (impl_->A_II_inv->apply(impl_->temp_I_local, impl_->x_I_local)) {
            // return false;
        }

        impl_->temp_I.restore_local_vector(impl_->temp_I_local);
        impl_->x_I.restore_local_vector(impl_->x_I_local);

        y = impl_->A_GG * impl_->x_G;
        impl_->x_G = impl_->A_GI * impl_->x_I;
        y -= impl_->x_G;
        return true;
    }

    bool SchurComplement<PetscMatrix>::finalize_from_eliminated_rhs(PetscVector &eliminated_rhs,
                                                                    const PetscVector &x_restricted,
                                                                    PetscVector &x) {
        UTOPIA_TRACE_REGION_BEGIN("SchurComplement::finalize");

        if (empty(x)) {
            UTOPIA_TRACE_REGION_END("SchurComplement::finalize");
            Utopia::Abort("SchurComplement<PetscMatrix>::finalize: x has to be initialized outside!");
        }

        impl_->x_I = impl_->A_IG * x_restricted;
        eliminated_rhs -= impl_->x_I;

        if (empty(impl_->x_I)) {
            impl_->x_I.zeros(layout(impl_->temp_I));
        } else {
            impl_->x_I.set(0);
        }

        eliminated_rhs.create_local_vector(impl_->temp_I_local);
        impl_->x_I.create_local_vector(impl_->x_I_local);

        assert(impl_->temp_I_local.comm().size() == 1);
        assert(impl_->x_I_local.comm().size() == 1);

        impl_->A_II_inv->apply(impl_->temp_I_local, impl_->x_I_local);

        eliminated_rhs.restore_local_vector(impl_->temp_I_local);
        impl_->x_I.restore_local_vector(impl_->x_I_local);

        impl_->restore_vector_blocks(x, impl_->x_I, x_restricted);

        UTOPIA_TRACE_REGION_END("SchurComplement::finalize");
        return true;
    }

    bool SchurComplement<PetscMatrix>::finalize(const PetscVector &rhs,
                                                const PetscVector &x_restricted,
                                                PetscVector &x) {
        UTOPIA_TRACE_REGION_BEGIN("SchurComplement::finalize");

        if (empty(x)) {
            UTOPIA_TRACE_REGION_END("SchurComplement::finalize");
            Utopia::Abort("SchurComplement<PetscMatrix>::finalize: x has to be initialized outside!");
        }

        bool ok = impl_->create_parallel_vector_block(rhs, impl_->is_eliminated, impl_->temp_I);

        finalize_from_eliminated_rhs(impl_->temp_I, x_restricted, x);

        UTOPIA_TRACE_REGION_END("SchurComplement::finalize");
        return ok;
    }

    Size SchurComplement<PetscMatrix>::size() const { return impl_->A_GG.size(); }

    Size SchurComplement<PetscMatrix>::local_size() const { return impl_->A_GG.local_size(); }

    SchurComplement<PetscMatrix>::Communicator &SchurComplement<PetscMatrix>::comm() { return impl_->A_GG.comm(); }

    const SchurComplement<PetscMatrix>::Communicator &SchurComplement<PetscMatrix>::comm() const {
        return impl_->A_GG.comm();
    }

    std::shared_ptr<PetscMatrix> SchurComplement<PetscMatrix>::reduced_matrix() const { return make_ref(impl_->A_GG); }

    void SchurComplement<PetscMatrix>::select(const PetscVector &x, PetscVector &x_free) const {
        impl_->create_parallel_vector_block(x, impl_->is_dof, x_free);
    }

}  // namespace utopia
