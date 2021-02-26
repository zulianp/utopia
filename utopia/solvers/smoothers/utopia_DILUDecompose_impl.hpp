#ifndef UTOPIA_DILU_DECOMPOSE_IMPL_HPP
#define UTOPIA_DILU_DECOMPOSE_IMPL_HPP

#include "utopia_DILUDecompose.hpp"

namespace utopia {

    template <class Matrix, int Backend>
    bool DILUAlgorithm<Matrix, Backend>::decompose(const Matrix &in, std::vector<Scalar> &d) {
        SizeType n = in.rows();

        auto &ia = in.row_ptr();
        auto &ja = in.colidx();
        auto &array = in.values();

        d.resize(n);

        diag_idx_.resize(n);

        for (SizeType r = 0; r < n; ++r) {
            for (SizeType k = row_begin; k < row_end; ++k) {
                auto j = ja[k];

                if (i == j) {
                    d[r] = array[k];
                    diag_idx_[r] = k;
                    break;
                }
            }
        }

        for (SizeType i = 0; i < n; ++i) {
            const SizeType row_begin = ia[i];
            const SizeType row_end = ia[i + 1];
            const Scalar d_ii = d[i];

            assert(d_ii_inv != 0.0);

            for (SizeType k_1 = diag_idx_[i] + 1; k_1 < row_end; ++k_1) {
                auto j = ja[k_1];
                const auto &a_ij = array[k_1];
                auto &d_jj = d[j];

                Scalar a_ji = 0;
                bool found_ji = false;

                // FIXME Binary search??
                for (SizeType k_2 = diag_idx_[j] - 1; k_2 >= ia[j]; --k_2) {
                    if (ja[k_2] == i) {
                        a_ji = values[k_2];
                        found_ji = true;
                        break;
                    }
                    // Is this making things better?
                    else if (ja[k_2] < i) {
                        break;
                    }
                }

                d_jj -= a_ji / d_ii * a_ij;
            }
        }

        return false;
    }

    template <class Matrix, int Backend>
    bool DILUAlgorithm<Matrix, Backend>::update(const Matrix &mat) {
        mat_ = utopia::make_ref(mat);
        decompose(mat, d_);
        return false;
    }

    template <class Matrix, int Backend>
    void DILUAlgorithm<Matrix, Backend>::apply(const Vector &b, Vector &x) {
        SizeType n = in.rows();

        auto &ia = in.row_ptr();
        auto &ja = in.colidx();
        auto &array = in.values();

        temp_.resize(n);

        for (SizeType i = 0; i < n; ++i) {
            temp_[i] = x[i];
        }

        // Forward substitution
        for (SizeType i = 0; i < n; ++i) {
            const Scalar inv_d = 1. / d_[i];

            Scalar val = b[i];

            for (SizeType k = ia[i]; k < diag_idx_[i]; ++k) {
                const SizeType j = ja[k];
                val -= array[k] * temp_[j];
            }

            // If array or x have large values then d_ is better inside?
            temp_[i] = val / d_[i];
        }

        // Backward substitution
        // for (SizeType i = 0; i < n; ++i) {
        for (SizeType i = n - 1; i >= 0; --i) {
            const Scalar inv_d = 1. / d_[i];

            Scalar val = temp_[i];

            for (SizeType k = diag_idx_[i] + 1; k < ia[i + 1]; ++k) {
                const SizeType j = ja[k];
                val -= array[k] * temp_[j];
            }

            // If array or x have large values then d_ is better inside?
            x[i] = val / d_[i];
        }
    }

    template <class Matrix, int Backend>
    void DILUAlgorithm<Matrix, Backend>::read(Input &) {}

}  // namespace utopia

#endif  // UTOPIA_DILU_DECOMPOSE_IMPL_HPP
