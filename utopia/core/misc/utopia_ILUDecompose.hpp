#ifndef UTOPIA_ILU_DECOMPOSE_HPP
#define UTOPIA_ILU_DECOMPOSE_HPP

#include "utopia_Algorithms.hpp"
#include "utopia_Views.hpp"

#include "utopia_CRSToBlockCRS.hpp"
#include "utopia_ILU.hpp"

#include <vector>

namespace utopia {

    template <typename SizeType>
    class DiagIdx {
    public:
        void init(SizeType n, const SizeType *ia, const SizeType *ja) {
            idx.resize(n);
            this->ia = ia;
            this->ja = ja;

            for (SizeType r = 0; r < n - 1; ++r) {
                const SizeType row_end = ia[r + 1];

                SizeType k;
                bool has_diag = false;
                for (k = ia[r]; k < row_end; ++k) {
                    if (ja[k] == r) {
                        has_diag = true;
                        break;
                    }
                }

                assert(has_diag);

                idx[r] = k;
            }
        }

        auto find_before_diag(const SizeType i, const SizeType j) const -> SizeType {
            SizeType search_end = idx[i];
            SizeType row_begin = ia[i];

            // maybe binary search?
            for (SizeType k = row_begin; k < search_end; ++k) {
                if (ja[k] == j) {
                    return k;
                }
            }

            return -1;
        };

        auto find_after_diag(const SizeType i, const SizeType j) const -> SizeType {
            SizeType search_begin = idx[i];
            SizeType row_end = ia[i + 1];

            // maybe binary search?
            for (SizeType k = search_begin; k < row_end; ++k) {
                if (ja[k] == j) {
                    return k;
                }
            }

            return -1;
        };

        auto find(const SizeType i, const SizeType j) const -> SizeType {
            if (i == j) return idx[i];

            if (j < i) {
                return find_before_diag(i, j);
            } else {
                return find_after_diag(i, j);
            }
        };

        std::vector<SizeType> idx;
        const SizeType *ia;
        const SizeType *ja;
    };

    template <class ScalarView, class IndexView>
    class ILUDecompose<CRSMatrix<ScalarView, IndexView, 1>, HOMEMADE> {
    public:
        using CRSMatrix = utopia::CRSMatrix<ScalarView, IndexView, 1>;
        using Scalar = typename std::remove_const<typename CRSMatrix::Scalar>::type;
        using SizeType = typename std::remove_const<typename CRSMatrix::SizeType>::type;
        using DiagIdx = utopia::DiagIdx<SizeType>;

        static bool decompose(CRSMatrix &in_out, const bool modified) {
            if (modified) {
                return decompose_aux(
                    [](const SizeType ii, const SizeType ij, const SizeType rj, const Scalar &e, Scalar *a) {
                        if (ij != -1) {
                            a[ij] -= e * a[rj];
                        } else {
                            a[ii] -= e * a[rj];
                        }
                    },
                    in_out);
            } else {
                return decompose_aux(
                    [](const SizeType, const SizeType ij, const SizeType rj, const Scalar &e, Scalar *a) {
                        if (ij != -1) {
                            a[ij] -= e * a[rj];
                        }
                    },
                    in_out);
            }
        }

        template <class ConstVector, class Vector>
        static void apply(const CRSMatrix &ilu, const ConstVector &b, Vector &L_inv_b, Vector &x) {
            device::fill(0.0, L_inv_b);
            device::fill(0.0, x);

            SizeType n = ilu.rows();
            auto &ia = ilu.row_ptr();
            auto &ja = ilu.colidx();
            auto &array = ilu.values();

            DiagIdx idx;
            idx.init(n, &ia[0], &ja[0]);

            // Forward substitution
            // https://algowiki-project.org/en/Forward_substitution#:~:text=Forward%20substitution%20is%20the%20process,math%5DL%5B%2Fmath%5D.
            for (SizeType i = 0; i < n; ++i) {
                const SizeType row_begin = ia[i];
                const SizeType row_diag = idx.idx[i];

                Scalar val = b[i];

                for (SizeType k = row_begin; k < row_diag; ++k) {
                    auto j = ja[k];
                    val -= array[k] * L_inv_b[j];
                }

                L_inv_b[i] = val / array[row_diag];
            }

            // Backward substitution
            // https://algowiki-project.org/en/Backward_substitution#:~:text=Backward%20substitution%20is%20a%20procedure,is%20a%20lower%20triangular%20matrix.
            for (SizeType i = n - 1; i >= 0; --i) {
                const SizeType row_end = ia[i + 1];
                const SizeType row_diag = idx.idx[i];

                Scalar val = L_inv_b[i];

                for (SizeType k = row_diag + 1; k < row_end; ++k) {
                    auto j = ja[k];
                    val -= array[k] * x[j];
                }

                x[i] = val / array[row_diag];
            }
        }

    private:
        template <typename Apply_IJ>
        static bool decompose_aux(Apply_IJ apply_ij, CRSMatrix &in_out) {
            SizeType n = in_out.rows();

            auto &ia = in_out.row_ptr();
            auto &ja = in_out.colidx();
            auto &array = in_out.values();

            DiagIdx idx;
            idx.init(n, &ia[0], &ja[0]);

            for (SizeType r = 0; r < n - 1; ++r) {
                const SizeType row_end = ia[r + 1];
                const SizeType k = idx.idx[r];
                const Scalar d = 1. / array[k];

                for (SizeType ki = k + 1; ki < row_end; ++ki) {
                    auto i = ja[ki];

                    auto ir = idx.find_before_diag(i, r);
                    if (ir == -1) continue;

                    auto ii = idx.idx[i];

                    Scalar e = (array[ir] *= d);

                    for (SizeType rj = k + 1; rj < row_end; ++rj) {
                        auto j = ja[rj];

                        auto ij = idx.find(i, j);
                        apply_ij(ii, ij, rj, e, &array[0]);
                    }
                }
            }

            // FIXME handle errors
            return true;
        }
    };

    template <class ScalarView, class IndexView, int BlockSize>
    class ILUDecompose<CRSMatrix<ScalarView, IndexView, BlockSize>, HOMEMADE> {
    public:
        using CRSMatrix = utopia::CRSMatrix<ScalarView, IndexView, BlockSize>;
        using Scalar = typename std::remove_const<typename CRSMatrix::Scalar>::type;
        using SizeType = typename std::remove_const<typename CRSMatrix::SizeType>::type;
        using DiagIdx = utopia::DiagIdx<SizeType>;
        using Block = utopia::StaticMatrix<Scalar, BlockSize, BlockSize>;
        using ArrayViewT = utopia::ArrayView<Scalar, DYNAMIC_SIZE, DYNAMIC_SIZE>;
        using BlockView = utopia::TensorView<ArrayViewT, 2>;

        using ConstArrayViewT = utopia::ArrayView<const Scalar, DYNAMIC_SIZE, DYNAMIC_SIZE>;
        using ConstBlockView = utopia::TensorView<ConstArrayViewT, 2>;

        static const int BlockSize2 = BlockSize * BlockSize;

        static bool decompose(CRSMatrix &in_out, const bool modified) {
            UTOPIA_UNUSED(modified);

            SizeType n_blocks = in_out.rows();

            auto &ia = in_out.row_ptr();
            auto &ja = in_out.colidx();
            // auto &array = in_out.values();

            DiagIdx idx;
            idx.init(n_blocks, &ia[0], &ja[0]);

            Block d_inv, e;
            BlockView d, air, arj, aij;

            d.raw_type().set_size(BlockSize, BlockSize);
            air.raw_type().set_size(BlockSize, BlockSize);
            arj.raw_type().set_size(BlockSize, BlockSize);
            aij.raw_type().set_size(BlockSize, BlockSize);

            for (SizeType r = 0; r < n_blocks - 1; ++r) {
                const SizeType row_end = ia[r + 1];
                const SizeType k = idx.idx[r];
                d.raw_type().set_data(in_out.block(k));
                assert(std::abs(det(d)) > 0);
                d_inv = inv(d);

                for (SizeType ki = k + 1; ki < row_end; ++ki) {
                    auto i = ja[ki];

                    auto ir = idx.find_before_diag(i, r);
                    if (ir == -1) continue;

                    air.raw_type().set_data(in_out.block(ir));

                    e = (d_inv * air);
                    air.copy(e);

                    for (SizeType rj = k + 1; rj < row_end; ++rj) {
                        auto j = ja[rj];

                        auto ij = idx.find(i, j);

                        arj.raw_type().set_data(in_out.block(rj));

                        if (ij != -1) {
                            aij.raw_type().set_data(in_out.block(ij));
                            aij -= e * arj;
                        }
                    }
                }
            }

            return true;
        }

        template <class ConstVector, class Vector>
        static void apply(const CRSMatrix &ilu, const ConstVector &b, Vector &L_inv_b, Vector &x) {
            device::fill(0.0, L_inv_b);
            device::fill(0.0, x);

            SizeType n_blocks = ilu.rows();
            auto &ia = ilu.row_ptr();
            auto &ja = ilu.colidx();

            DiagIdx idx;
            idx.init(n_blocks, &ia[0], &ja[0]);

            Block d_inv;
            ConstBlockView d;

            d.raw_type().set_size(BlockSize, BlockSize);

            StaticVector<Scalar, BlockSize> val;

            // Forward substitution
            for (SizeType i = 0; i < n_blocks; ++i) {
                const SizeType i_offset = i * BlockSize;
                const SizeType row_begin = ia[i];
                const SizeType row_diag = idx.idx[i];

                d.raw_type().set_data(ilu.block(row_diag));
                d_inv = inv(d);

                for (SizeType sub_i = 0; sub_i < BlockSize; ++sub_i) {
                    val[sub_i] = b[i_offset + sub_i];
                }

                for (SizeType k = row_begin; k < row_diag; ++k) {
                    const SizeType j = ja[k];
                    const SizeType j_offset = j * BlockSize;

                    auto *block = ilu.block(k);

                    // Matrix-vector multiplication (in the block)
                    for (SizeType sub_i = 0; sub_i < BlockSize; ++sub_i) {
                        const SizeType k_offset_sub_i = sub_i * BlockSize;

                        for (SizeType bj = 0; bj < BlockSize; ++bj) {
                            val[sub_i] -= block[k_offset_sub_i + bj] * L_inv_b[j_offset + bj];
                        }
                    }
                }

                auto expr = d_inv * val;
                for (SizeType sub_i = 0; sub_i < BlockSize; ++sub_i) {
                    L_inv_b[i_offset + sub_i] = expr(sub_i);
                }
            }

            // // Backward substitution
            for (SizeType i = n_blocks - 1; i >= 0; --i) {
                const SizeType i_offset = i * BlockSize;
                const SizeType row_diag = idx.idx[i];
                const SizeType row_end = ia[i + 1];

                d.raw_type().set_data(ilu.block(row_diag));
                d_inv = inv(d);

                for (SizeType sub_i = 0; sub_i < BlockSize; ++sub_i) {
                    val[sub_i] = L_inv_b[i_offset + sub_i];
                }

                for (SizeType k = row_diag + 1; k < row_end; ++k) {
                    const SizeType j = ja[k];
                    const SizeType j_offset = j * BlockSize;

                    auto *block = ilu.block(k);

                    // Matrix-vector multiplication (in the block)
                    for (SizeType sub_i = 0; sub_i < BlockSize; ++sub_i) {
                        const SizeType k_offset_sub_i = sub_i * BlockSize;

                        for (SizeType sub_j = 0; sub_j < BlockSize; ++sub_j) {
                            val[sub_i] -= block[k_offset_sub_i + sub_j] * x[j_offset + sub_j];
                        }
                    }
                }

                auto expr = d_inv * val;
                for (SizeType sub_i = 0; sub_i < BlockSize; ++sub_i) {
                    x[i_offset + sub_i] = expr(sub_i);
                }
            }
        }
    };  // namespace utopia

    template <class Matrix>
    bool ilu_decompose(Matrix &mat, const bool modified = false) {
        return ILUDecompose<Matrix>::decompose(mat, modified);
    }

    template <class Matrix, class ConstVector, class Vector>
    void ilu_apply(const Matrix &ilu, const ConstVector &b, Vector &L_inv_b, Vector &out) {
        ILUDecompose<Matrix>::apply(ilu, b, L_inv_b, out);
    }

}  // namespace utopia

#endif  // UTOPIA_ILU_DECOMPOSE_HPP