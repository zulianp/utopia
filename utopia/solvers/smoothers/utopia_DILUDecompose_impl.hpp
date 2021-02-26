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
                auto j = ja[k_1];
                const auto &a_ij = array[k_1];
                auto &d_jj = d[j];

                auto k_2 = transpose_idx_.transpose(k_1);

                if (k_2 != -1) {
                    d_jj -= (array[k_2] / d_ii) * a_ij;
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

        for (SizeType i = 0; i < n; ++i) {
            L_inv_b_[i] = x[i];
            // x[i] = 0.0;
        }

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
                val -= array[k] * L_inv_b_[j];
                // val -= array[k] * x[j];
            }

            // If array or x have large values then d_ is better inside?
            x[i] = L_inv_b_[i] - val / d_[i];
        }

        // SizeType n = mat_->rows();

        // L_inv_b_.resize(n);
        // device::fill(0.0, L_inv_b_);
        // device::fill(0.0, x);

        // auto &&ia = mat_->row_ptr();
        // auto &&ja = mat_->colidx();
        // auto &&array = mat_->values();

        // // Forward substitution
        // //
        // https://algowiki-project.org/en/Forward_substitution#:~:text=Forward%20substitution%20is%20the%20process,math%5DL%5B%2Fmath%5D.
        // for (SizeType i = 0; i < n; ++i) {
        //     const SizeType row_begin = ia[i];
        //     const SizeType row_diag = diag_idx_.diag(i);

        //     assert(array[row_diag] != 0.0);

        //     Scalar val = b[i];
        //     assert(std::isfinite(val));

        //     for (SizeType k = row_begin; k < row_diag; ++k) {
        //         auto j = ja[k];
        //         assert(j < n);

        //         val -= array[k] * L_inv_b_[j];
        //     }

        //     L_inv_b_[i] = val / d_[i];

        //     assert(L_inv_b_[i] == L_inv_b_[i]);
        //     assert(std::isfinite(L_inv_b_[i]));
        // }

        // // Backward substitution
        // //
        // https://algowiki-project.org/en/Backward_substitution#:~:text=Backward%20substitution%20is%20a%20procedure,is%20a%20lower%20triangular%20matrix.
        // for (SizeType i = n - 1; i >= 0; --i) {
        //     const SizeType row_end = ia[i + 1];
        //     const SizeType row_diag = diag_idx_.diag(i);

        //     assert(array[row_diag] != 0.0);

        //     Scalar val = L_inv_b_[i];

        //     for (SizeType k = row_diag + 1; k < row_end; ++k) {
        //         auto j = ja[k];
        //         val -= array[k] * x[j];
        //     }

        //     x[i] = val / d_[i];
        //     assert(x[i] == x[i]);
        // }
    }

    template <class Matrix, class Vector>
    void DILUAlgorithm<Matrix, Vector>::read(Input &) {}

}  // namespace utopia

#endif  // UTOPIA_DILU_DECOMPOSE_IMPL_HPP
