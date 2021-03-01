#ifndef UTOPIA_DILU_DECOMPOSE_IMPL_HPP
#define UTOPIA_DILU_DECOMPOSE_IMPL_HPP

#include "utopia_DILUDecompose.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    bool DILUAlgorithm<Matrix, Vector>::update(const Matrix &mat) {
        UTOPIA_TRACE_REGION_BEGIN("DILUAlgorithm::update");

        mat_ = utopia::make_ref(mat);

        SizeType n = mat_->rows();

        auto &&ia = mat_->row_ptr();
        auto &&ja = mat_->colidx();
        auto &&array = mat_->values();

        d_.resize(n);
        diag_idx_.init(n, &ia[0], &ja[0]);
        transpose_idx_.init(n, array.size(), &ia[0], &ja[0]);

        for (SizeType r = 0; r < n; ++r) {
            d_[r] = array[diag_idx_.diag(r)];
        }

        for (SizeType i = 0; i < n; ++i) {
            const SizeType row_end = ia[i + 1];
            const Scalar d_ii = d_[i];

            assert(d_ii != 0.0);

            for (SizeType k_1 = diag_idx_.diag(i) + 1; k_1 < row_end; ++k_1) {
                const SizeType j = ja[k_1];
                const Scalar &a_ij = array[k_1];
                Scalar &d_jj = d_[j];

                SizeType k_2 = transpose_idx_.transpose(k_1);

                if (k_2 != -1) {
                    const Scalar a_ji = array[k_2];
                    d_jj -= (a_ji / d_ii) * a_ij;
                }
            }
        }

        UTOPIA_TRACE_REGION_END("DILUAlgorithm::update");
        return false;
    }

    // Works with laplacian
    template <class Matrix, class Vector>
    void DILUAlgorithm<Matrix, Vector>::apply(const Vector &b, Vector &x) {
        UTOPIA_TRACE_REGION_BEGIN("DILUAlgorithm::apply");

        const SizeType n = mat_->rows();

        auto &&ia = mat_->row_ptr();
        auto &&ja = mat_->colidx();
        auto &&array = mat_->values();

        L_inv_b_.resize(n);

        device::fill(0.0, L_inv_b_);

        // Forward substitution
        for (SizeType i = 0; i < n; ++i) {
            Scalar val = b[i];

            for (SizeType k = ia[i]; k < diag_idx_.diag(i); ++k) {
                const SizeType j = ja[k];
                val -= array[k] * L_inv_b_[j];
            }

            // If array or x have large values then d_ is better inside?
            L_inv_b_[i] = val / d_[i];
        }

        // Backward substitution
        for (SizeType i = n - 1; i >= 0; --i) {
            Scalar val = 0.0;

            for (SizeType k = diag_idx_.diag(i) + 1; k < ia[i + 1]; ++k) {
                const SizeType j = ja[k];
                val += array[k] * x[j];
            }

            // If array or x have large values then d_ is better inside?
            x[i] = L_inv_b_[i] - val / d_[i];
        }

        UTOPIA_TRACE_REGION_END("DILUAlgorithm::apply");
    }

    template <class Matrix, class Vector>
    void DILUAlgorithm<Matrix, Vector>::read(Input &) {}

    //////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    template <class Matrix, class Vector, int BlockSize>
    bool BlockDILUAlgorithm<Matrix, Vector, BlockSize>::update(const Matrix &mat) {
        UTOPIA_TRACE_REGION_BEGIN("BlockDILUAlgorithm::update");

        mat_ = utopia::make_ref(mat);

        SizeType n_blocks = mat_->rows();

        auto &&ia = mat_->row_ptr();
        auto &&ja = mat_->colidx();
        // auto &&array = mat_->values();

        d_.resize(n_blocks);

        diag_idx_.init(n_blocks, &ia[0], &ja[0]);
        transpose_idx_.init(n_blocks, mat_->n_blocks(), &ia[0], &ja[0]);

        BlockView d_view, a_ij, a_ji;
        d_view.raw_type().set_size(BlockSize, BlockSize);
        a_ij.raw_type().set_size(BlockSize, BlockSize);
        a_ji.raw_type().set_size(BlockSize, BlockSize);

        Block d_ii_inv;

        for (SizeType i = 0; i < n_blocks; ++i) {
            d_view.raw_type().set_data(mat_->block(i));
            d_[i].copy(d_view);
        }

        for (SizeType block_i = 0; block_i < n_blocks; ++block_i) {
            const SizeType row_end = ia[block_i + 1];

            const auto &d_ii = d_[block_i];
            d_ii_inv = inv(d_ii);

            for (SizeType k_1 = diag_idx_.diag(block_i) + 1; k_1 < row_end; ++k_1) {
                const SizeType block_j = ja[k_1];
                a_ij.raw_type().set_data(mat_->block(k_1));

                auto &d_jj = d_[block_j];

                SizeType k_2 = transpose_idx_.transpose(k_1);

                if (k_2 != -1) {
                    a_ji.raw_type().set_data(mat_->block(k_2));
                    d_jj -= a_ji * d_ii_inv * a_ij;
                }
            }
        }

        Block temp;
        for (SizeType i = 0; i < n_blocks; ++i) {
            temp.copy(d_[i]);
            d_[i] = inv(temp);
        }

        UTOPIA_TRACE_REGION_END("BlockDILUAlgorithm::update");
        return false;
    }

    // Works with laplacian
    template <class Matrix, class Vector, int BlockSize>
    void BlockDILUAlgorithm<Matrix, Vector, BlockSize>::apply(const Vector &b, Vector &x) {
        UTOPIA_TRACE_REGION_BEGIN("BlockDILUAlgorithm::apply");

        const SizeType n_blocks = mat_->rows();

        auto &&ia = mat_->row_ptr();
        auto &&ja = mat_->colidx();
        // auto &&array = mat_->values();

        L_inv_b_.resize(n_blocks);

        device::fill(0.0, L_inv_b_);

        StaticVector<Scalar, BlockSize> val, x_vec;

        BlockView a_ij;
        a_ij.raw_type().set_size(BlockSize, BlockSize);

        // Forward substitution
        for (SizeType block_i = 0; block_i < n_blocks; ++block_i) {
            for (SizeType sub_i = 0; sub_i < BlockSize; ++sub_i) {
                val[sub_i] = b[block_i * BlockSize + sub_i];
            }

            for (SizeType block_k = ia[block_i]; block_k < diag_idx_.diag(block_i); ++block_k) {
                const SizeType block_j = ja[block_k];

                a_ij.raw_type().set_data(mat_->block(block_k));

                for (int sub_j = 0; sub_j < BlockSize; ++sub_j) {
                    x_vec[sub_j] = L_inv_b_[block_j * BlockSize + sub_j];
                }

                val -= a_ij * x_vec;
            }

            auto expr = d_[block_i] * val;
            for (SizeType sub_i = 0; sub_i < BlockSize; ++sub_i) {
                L_inv_b_[block_i * BlockSize + sub_i] = expr(sub_i);
            }
        }

        // Backward substitution
        for (SizeType block_i = n_blocks - 1; block_i >= 0; --block_i) {
            for (SizeType sub_i = 0; sub_i < BlockSize; ++sub_i) {
                val[sub_i] = 0.0;
            }

            for (SizeType block_k = diag_idx_.diag(block_i) + 1; block_k < ia[block_i + 1]; ++block_k) {
                const SizeType block_j = ja[block_k];

                a_ij.raw_type().set_data(mat_->block(block_k));

                for (int sub_j = 0; sub_j < BlockSize; ++sub_j) {
                    x_vec[sub_j] = x[block_j * BlockSize + sub_j];
                }

                val += a_ij * x_vec;
            }

            // If array or x have large values then d_ is better inside?
            auto expr = d_[block_i] * val;
            for (SizeType sub_i = 0; sub_i < BlockSize; ++sub_i) {
                const SizeType i = block_i * BlockSize + sub_i;
                x[i] = L_inv_b_[i] - expr(sub_i);
            }
        }

        UTOPIA_TRACE_REGION_END("BlockDILUAlgorithm::apply");
    }

    template <class Matrix, class Vector, int BlockSize>
    void BlockDILUAlgorithm<Matrix, Vector, BlockSize>::read(Input &) {}

}  // namespace utopia

#endif  // UTOPIA_DILU_DECOMPOSE_IMPL_HPP
