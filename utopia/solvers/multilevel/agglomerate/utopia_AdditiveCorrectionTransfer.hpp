#ifndef UTOPIA_ADDITIVE_CORRECTION_TRANSFER_HPP
#define UTOPIA_ADDITIVE_CORRECTION_TRANSFER_HPP

#include "utopia_Transfer.hpp"

namespace utopia {

    template <class Matrix, class Vector, int BlockSize = 1>
    class AdditiveCorrectionTransfer final : public Transfer<Matrix, Vector> {
    public:
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        using IndexArray = typename utopia::Traits<Vector>::IndexArray;

        bool interpolate(const Vector &x_coarse, Vector &x_fine) const override {
            UTOPIA_UNUSED(x_coarse);
            UTOPIA_UNUSED(x_fine);
            // auto vl = layout(x_coarse);
            // assert(vl.local_size(0) == n_coarse_local_);

            // assert(!empty(x_fine));

            // auto x_coarse_view = local_view_device(x_coarse);
            // auto x_fine_view = local_view_device(x_fine);

            // for (SizeType i = 0; i < n_coarse_local_; ++i) {
            //     x_fine_view.set(i, x_coarse_view.get(parent_[i]));
            // }

            return false;
        }

        bool restrict(const Vector &x_fine, Vector &x_coarse) const override {
            UTOPIA_UNUSED(x_coarse);
            UTOPIA_UNUSED(x_fine);
            // auto vl = layout(x_fine);

            // const SizeType n_fine = vl.size(0);

            // if (empty(x_coarse)) {
            //     x_coarse.zeros(layout(x_fine.comm(), n_coarse_local_, n_coarse_global_));
            // } else {
            //     x_coarse.set(0.0);
            // }

            // auto x_fine_view = local_view_device(x_fine);
            // auto x_coarse_view = local_view_device(x_coarse);
            // for (SizeType i = 0; i < n_fine; ++i) {
            //     x_coarse_view.add(parent_[i], x_fine_view.get(i));
            // }

            // return true;

            return false;
        }

        bool restrict(const Matrix &mat_fine, Matrix &mat_coarse) const override {
            UTOPIA_UNUSED(mat_coarse);
            UTOPIA_UNUSED(mat_fine);
            assert(false);

            // IndexArray d_nnz(n_coarse_local_, 0);
            // // // FIXME when doing a proper parallel version
            // IndexArray o_nnz(n_coarse_local_, 0);

            // auto crs_fine = crs_view(mat_fine);

            // std::vector<bool> touched(n_coarse_local_, false);
            // const SizeType n_fine = mat_fine.local_rows();

            // for (SizeType i = 0; i < n_fine; ++i) {
            //     const auto &&row = crs_fine.row(i);
            //     const SizeType n_values = row.n_blocks();
            //     const SizeType parent_i = parent_[i];

            //     for (SizeType k = 0; k < n_values; ++k) {
            //         const SizeType j = row.colidx(k);
            //         const SizeType parent_j = parent_[j];

            //         if (!touched[parent_j]) {
            //             touched[parent_j] = true;
            //             ++d_nnz[parent_i];
            //         }
            //     }

            //     for (SizeType k = 0; k < n_values; ++k) {
            //         const SizeType j = row.colidx(k);
            //         const SizeType parent_j = parent_[j];
            //         touched[parent_j] = false;
            //     }
            // }

            // if (BlockSize == 1) {
            //     mat_coarse.sparse(
            //         layout(mat_fine.comm(), n_coarse_local_, n_coarse_local_, n_coarse_global_, n_coarse_global_),
            //         d_nnz,
            //         o_nnz);

            // } else {
            //     assert(false && "IMPLEMENT ME");
            //     mat_coarse.block_sparse(
            //         layout(mat_fine.comm(), n_coarse_local_, n_coarse_local_, n_coarse_global_, n_coarse_global_),
            //         d_nnz,
            //         o_nnz,
            //         BlockSize);
            // }

            // Write<Matrix> w(mat_coarse);
            // for (SizeType fine_i = 0; fine_i < n_fine; ++fine_i) {
            //     auto row = crs_fine.row(fine_i);
            //     const SizeType n_values = row.n_blocks();

            //     const SizeType coarse_i = parent_[fine_i];

            //     for (SizeType k = 0; k < n_values; ++k) {
            //         const SizeType fine_j = row.colidx(k);
            //         const SizeType coarse_j = parent_[fine_j];

            //         const Scalar val = row.value(k);
            //         // add to offdiagonals
            //         mat_coarse.c_add(coarse_i, coarse_j, val);

            //         // add also to diagonal
            //         mat_coarse.c_add(coarse_i, coarse_i, val);
            //     }
            // }

            return false;
        }

        bool boolean_restrict_or(const Vector &, Vector &) override {
            assert(false && "IMPLENT ME");
            return false;
        }

        bool project_down(const Vector &, Vector &) const override {
            assert(false && "IMPLENT ME");
            return false;
        }

        bool project_down_positive_negative(const Vector &, const Vector &, Vector &) override {
            assert(false && "IMPLENT ME");
            return false;
        }

        void init_memory() override {}
        Scalar interpolation_inf_norm() const override {
            assert(false && "IMPLENT ME");
            return -1.0;
        }
        Scalar projection_inf_norm() const override {
            assert(false && "IMPLENT ME");
            return -1.0;
        }
        Scalar restriction_inf_norm() const override {
            assert(false && "IMPLENT ME");
            return -1.0;
        }

        void handle_equality_constraints(const Vector &) override { assert(false && "IMPLENT ME"); }

        IndexArray &parent() { return parent_; }

        void set_size(const SizeType n_coarse_local, SizeType n_coarse_global) {
            n_coarse_local_ = n_coarse_local;
            n_coarse_global_ = n_coarse_global;
        }

    private:
        SizeType n_coarse_local_{-1}, n_coarse_global_{-1};
        IndexArray parent_;
    };
}  // namespace utopia

#endif  // UTOPIA_ADDITIVE_CORRECTION_TRANSFER_HPP
