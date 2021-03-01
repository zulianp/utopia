#ifndef UTOPIA_DILU_DECOMPOSE_IMPL_HPP
#define UTOPIA_DILU_DECOMPOSE_IMPL_HPP

#include "utopia_DILUDecompose.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    bool DILUAlgorithm<Matrix, Vector>::decompose(const Matrix &in, std::vector<Scalar> &d) {
        SizeType n = in.rows();

        auto &&ia = in.row_ptr();
        auto &&ja = in.colidx();
        auto &&array = in.values();

        d.resize(n);
        diag_idx_.init(n, &ia[0], &ja[0]);
        transpose_idx_.init(n, array.size(), &ia[0], &ja[0]);

        for (SizeType r = 0; r < n; ++r) {
            d[r] = array[diag_idx_.diag(r)];
        }

        for (SizeType i = 0; i < n; ++i) {
            const SizeType row_end = ia[i + 1];
            const Scalar d_ii = d[i];

            assert(d_ii != 0.0);

            for (SizeType k_1 = diag_idx_.diag(i) + 1; k_1 < row_end; ++k_1) {
                const SizeType j = ja[k_1];
                const Scalar &a_ij = array[k_1];
                Scalar &d_jj = d[j];

                SizeType k_2 = transpose_idx_.transpose(k_1);

                if (k_2 != -1) {
                    const Scalar a_ji = array[k_2];
                    d_jj -= (a_ji / d_ii) * a_ij;
                }
            }
        }

        return false;
    }

    template <class Matrix, class Vector>
    bool DILUAlgorithm<Matrix, Vector>::update(const Matrix &mat) {
        mat_ = utopia::make_ref(mat);
        decompose(mat, d_);
        return false;
    }

    template <class Matrix, class Vector>
    void DILUAlgorithm<Matrix, Vector>::apply(const Vector &b, Vector &x) {
        const SizeType n = mat_->rows();

        auto &&ia = mat_->row_ptr();
        auto &&ja = mat_->colidx();
        auto &&array = mat_->values();

        L_inv_b_.resize(n);

        device::fill(0.0, L_inv_b_);
        device::fill(0.0, x);

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
    }

    template <class Matrix, class Vector>
    void DILUAlgorithm<Matrix, Vector>::read(Input &) {}

}  // namespace utopia

#endif  // UTOPIA_DILU_DECOMPOSE_IMPL_HPP
