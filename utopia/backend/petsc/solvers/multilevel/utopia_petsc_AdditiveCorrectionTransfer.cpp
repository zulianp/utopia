#include "utopia_petsc_AdditiveCorrectionTransfer.hpp"

#include "utopia_petsc_CrsView.hpp"

#include "utopia_DeviceView.hpp"
#include "utopia_Writable.hpp"

namespace utopia {

    bool AdditiveCorrectionTransfer<PetscMatrix, PetscVector>::interpolate(const PetscVector &x_coarse,
                                                                           PetscVector &x_fine) const {
        auto vl = layout(x_coarse);
        assert(vl.local_size(0) == n_coarse_local_);

        assert(!empty(x_fine));

        auto x_coarse_view = local_view_device(x_coarse);
        auto x_fine_view = local_view_device(x_fine);

        for (SizeType i = 0; i < n_coarse_local_; ++i) {
            x_fine_view.set(i, x_coarse_view.get(parent_[i]));
        }

        return false;
    }

    bool AdditiveCorrectionTransfer<PetscMatrix, PetscVector>::restrict(const PetscVector &x_fine,
                                                                        PetscVector &x_coarse) const {
        auto vl = layout(x_fine);

        const SizeType n_fine = vl.size(0);

        if (empty(x_coarse)) {
            x_coarse.zeros(layout(x_fine.comm(), n_coarse_local_, n_coarse_global_));
        } else {
            x_coarse.set(0.0);
        }

        auto x_fine_view = local_view_device(x_fine);
        auto x_coarse_view = local_view_device(x_coarse);
        for (SizeType i = 0; i < n_fine; ++i) {
            x_coarse_view.add(parent_[i], x_fine_view.get(i));
        }

        return true;
    }

    bool AdditiveCorrectionTransfer<PetscMatrix, PetscVector>::restrict(const PetscMatrix &mat_fine,
                                                                        PetscMatrix &mat_coarse) const {
        IndexArray d_nnz(n_coarse_local_, 0);
        // // FIXME when doing a proper parallel version
        IndexArray o_nnz(n_coarse_local_, 0);

        auto crs_fine = crs_view(mat_fine);

        std::vector<bool> touched(n_coarse_local_, false);
        const SizeType n_fine = mat_fine.local_rows();

        for (SizeType i = 0; i < n_fine; ++i) {
            const auto &&row = crs_fine.row(i);
            const SizeType n_values = row.length;
            const SizeType parent_i = parent_[i];

            for (SizeType k = 0; k < n_values; ++k) {
                const SizeType j = row.colidx(k);
                const SizeType parent_j = parent_[j];

                if (!touched[parent_j]) {
                    touched[parent_j] = true;
                    ++d_nnz[parent_i];
                }
            }

            for (SizeType k = 0; k < n_values; ++k) {
                const SizeType j = row.colidx(k);
                const SizeType parent_j = parent_[j];
                touched[parent_j] = false;
            }
        }

        mat_coarse.sparse(layout(mat_fine.comm(), n_coarse_local_, n_coarse_local_, n_coarse_global_, n_coarse_global_),
                          d_nnz,
                          o_nnz);

        Write<PetscMatrix> w(mat_coarse);
        for (SizeType fine_i = 0; fine_i < n_fine; ++fine_i) {
            auto row = crs_fine.row(fine_i);
            const SizeType n_values = row.length;

            const SizeType coarse_i = parent_[fine_i];

            for (SizeType k = 0; k < n_values; ++k) {
                const SizeType fine_j = row.colidx(k);
                const SizeType coarse_j = parent_[fine_j];

                const PetscScalar val = row.value(k);
                // add to offdiagonals
                mat_coarse.c_add(coarse_i, coarse_j, val);

                // add also to diagonal
                mat_coarse.c_add(coarse_i, coarse_i, val);
            }
        }

        return false;
    }

    bool AdditiveCorrectionTransfer<PetscMatrix, PetscVector>::boolean_restrict_or(const PetscVector &, PetscVector &) {
        assert(false && "IMPLENT ME");
        return false;
    }

    bool AdditiveCorrectionTransfer<PetscMatrix, PetscVector>::project_down(const PetscVector &, PetscVector &) const {
        assert(false && "IMPLENT ME");
        return false;
    }

    bool AdditiveCorrectionTransfer<PetscMatrix, PetscVector>::project_down_positive_negative(const PetscVector &,
                                                                                              const PetscVector &,
                                                                                              PetscVector &) {
        assert(false && "IMPLENT ME");
        return false;
    }

    void AdditiveCorrectionTransfer<PetscMatrix, PetscVector>::init_memory() {}

    PetscScalar AdditiveCorrectionTransfer<PetscMatrix, PetscVector>::interpolation_inf_norm() const {
        assert(false && "IMPLENT ME");
        return -1.0;
    }
    PetscScalar AdditiveCorrectionTransfer<PetscMatrix, PetscVector>::projection_inf_norm() const {
        assert(false && "IMPLENT ME");
        return -1.0;
    }
    PetscScalar AdditiveCorrectionTransfer<PetscMatrix, PetscVector>::restriction_inf_norm() const {
        assert(false && "IMPLENT ME");
        return -1.0;
    }

    void AdditiveCorrectionTransfer<PetscMatrix, PetscVector>::handle_equality_constraints(const PetscVector &) {
        assert(false && "IMPLENT ME");
    }

}  // namespace utopia
