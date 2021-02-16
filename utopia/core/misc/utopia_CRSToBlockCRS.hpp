#ifndef UTOPIA_CRS_TO_BLOCK_CRS_HPP
#define UTOPIA_CRS_TO_BLOCK_CRS_HPP

#include "utopia_Base.hpp"

#include "utopia_Algorithms.hpp"
#include "utopia_Traits.hpp"

#include "utopia_CRSMatrix.hpp"

#include <vector>

namespace utopia {

    template <class ScalarView,
              class IndexView,
              int BlockSize,
              class OutScalarView = ScalarView,
              class OutIndexView = IndexView>
    class CRSToBlockCRS {
    public:
        using Scalar = typename std::remove_const<typename Traits<ScalarView>::ValueType>::type;
        using SizeType = typename std::remove_const<typename Traits<IndexView>::ValueType>::type;

        static void apply(const CRSMatrix<ScalarView, IndexView, 1> &in,
                          CRSMatrix<OutScalarView, OutIndexView, BlockSize> &out,
                          const bool sort_columns = true) {
            apply(in.rows(), in.row_ptr(), in.colidx(), in.values(), out, sort_columns);
        }

        template <typename IA, typename JA, typename Array>
        static void apply(const SizeType n,
                          const IA &ia,
                          const JA &ja,
                          const Array &array,
                          CRSMatrix<OutScalarView, OutIndexView, BlockSize> &out,
                          const bool sort_columns = true) {
            static const int BlockSize_2 = BlockSize * BlockSize;
            auto n_blocks = n / BlockSize;
            assert(n_blocks * BlockSize == n);
            assert(n > 0);

            out.set_cols(n_blocks);

            std::vector<SizeType> block_pattern(n_blocks);

            auto &row_ptr = out.row_ptr();
            auto &colidx = out.colidx();
            auto &values = out.values();

            row_ptr.resize(n_blocks + 1, 0);

            for (SizeType block_i = 0; block_i < n_blocks; block_i++) {
                const SizeType offset_i = block_i * BlockSize;

                for (SizeType i = offset_i; i < offset_i + BlockSize; ++i) {
                    const SizeType row_begin = ia[i];
                    const SizeType row_end = ia[i + 1];

                    for (SizeType k = row_begin; k < row_end; ++k) {
                        SizeType j = ja[k];
                        SizeType block_j = j / BlockSize;

                        assert(block_j < SizeType(block_pattern.size()));

                        if (block_pattern[block_j] == 0) {
                            assert((block_i + 1) < SizeType(row_ptr.size()));
                            row_ptr[block_i + 1]++;
                            block_pattern[block_j] = 1.0;
                        }
                    }
                }

                for (SizeType i = offset_i; i < offset_i + BlockSize; ++i) {
                    const SizeType row_begin = ia[i];
                    const SizeType row_end = ia[i + 1];

                    for (SizeType k = row_begin; k < row_end; ++k) {
                        SizeType j = ja[k];
                        SizeType block_j = j / BlockSize;
                        block_pattern[block_j] = 0;
                    }
                }
            }

            for (SizeType i = 0; i < n_blocks; ++i) {
                row_ptr[i + 1] += row_ptr[i];
            }

            colidx.resize(row_ptr[n_blocks], 0);
            std::vector<SizeType> count(n_blocks, 0);

            for (SizeType block_i = 0; block_i < n_blocks; block_i++) {
                const SizeType offset_i = block_i * BlockSize;

                for (SizeType i = offset_i; i < offset_i + BlockSize; ++i) {
                    const SizeType row_begin = ia[i];
                    const SizeType row_end = ia[i + 1];

                    for (SizeType k = row_begin; k < row_end; ++k) {
                        SizeType j = ja[k];
                        SizeType block_j = j / BlockSize;

                        if (block_pattern[block_j] == 0) {
                            block_pattern[block_j] = 1.0;
                            SizeType idx_offset = row_ptr[block_i] + count[block_i];
                            colidx[idx_offset] = block_j;
                            ++count[block_i];
                        }
                    }
                }

                for (SizeType i = offset_i; i < offset_i + BlockSize; ++i) {
                    const SizeType row_begin = ia[i];
                    const SizeType row_end = ia[i + 1];

                    for (SizeType k = row_begin; k < row_end; ++k) {
                        SizeType j = ja[k];
                        SizeType block_j = j / BlockSize;
                        block_pattern[block_j] = 0;
                    }
                }

                if (sort_columns) {
                    std::sort(&colidx[row_ptr[block_i]], &colidx[row_ptr[block_i + 1] - 1]);
                }
            }

            values.resize(row_ptr[n_blocks] * BlockSize_2);

            update(n, ia, ja, array, out);
        }

        template <typename IA, typename JA, typename Array>
        static void update(const SizeType n,
                           const IA &ia,
                           const JA &ja,
                           const Array &array,
                           CRSMatrix<OutScalarView, OutIndexView, BlockSize> &out) {
            static const int BlockSize_2 = BlockSize * BlockSize;

            UTOPIA_UNUSED(n);

            const SizeType n_blocks = out.rows();
            assert(n == (BlockSize * n_blocks));

            auto &row_ptr = out.row_ptr();
            auto &colidx = out.colidx();
            auto &values = out.values();
            device::fill(0.0, values);

            for (SizeType block_i = 0; block_i < n_blocks; block_i++) {
                const SizeType offset_i = block_i * BlockSize;
                const SizeType block_row_begin = row_ptr[block_i];
                const SizeType block_row_end = row_ptr[block_i + 1];

                for (SizeType i = offset_i; i < offset_i + BlockSize; ++i) {
                    const SizeType row_begin = ia[i];
                    const SizeType row_end = ia[i + 1];
                    const SizeType sub_i = i - offset_i;

                    SizeType current_block = block_row_begin;

                    for (SizeType k = row_begin; k < row_end; ++k) {
                        const SizeType j = ja[k];
                        const Scalar value = array[k];

                        const SizeType block_j = j / BlockSize;
                        const SizeType sub_j = j - block_j * BlockSize;

                        SizeType k_block = current_block;
                        for (; k_block < block_row_end; ++k_block) {
                            if (colidx[k_block] == block_j) {
                                current_block = k_block;
                                break;
                            }
                        }

                        assert(colidx[k_block] == block_j);

                        values[k_block * BlockSize_2 + (sub_i * BlockSize) + sub_j] = value;
                    }
                }
            }
        }

        //////////////////////////////////////////////////////////////////////////////////////////

        template <typename DiagView>
        static void apply_split_diag(const CRSMatrix<ScalarView, IndexView, 1> &in,
                                     CRSMatrix<OutScalarView, OutIndexView, BlockSize> &out,
                                     DiagView &diag,
                                     const bool sort_columns = true) {
            apply_split_diag(in.rows(), in.row_ptr(), in.colidx(), in.values(), out, diag, sort_columns);
        }

        template <typename IA, typename JA, typename Array, typename Diag>
        static void apply_split_diag(const SizeType n,
                                     const IA &ia,
                                     const JA &ja,
                                     const Array &array,
                                     CRSMatrix<OutScalarView, OutIndexView, BlockSize> &out,
                                     Diag &diag,
                                     const bool sort_columns = true) {
            static const int BlockSize_2 = BlockSize * BlockSize;
            auto n_blocks = n / BlockSize;
            assert(n_blocks * BlockSize == n);
            assert(n > 0);

            out.set_cols(n_blocks);

            std::vector<SizeType> block_pattern(n_blocks);

            auto &row_ptr = out.row_ptr();
            auto &colidx = out.colidx();
            auto &values = out.values();

            row_ptr.resize(n_blocks + 1, 0);
            diag.resize(n_blocks * BlockSize_2);

            for (SizeType block_i = 0; block_i < n_blocks; block_i++) {
                const SizeType offset_i = block_i * BlockSize;

                for (SizeType i = offset_i; i < offset_i + BlockSize; ++i) {
                    const SizeType row_begin = ia[i];
                    const SizeType row_end = ia[i + 1];

                    for (SizeType k = row_begin; k < row_end; ++k) {
                        SizeType j = ja[k];
                        SizeType block_j = j / BlockSize;

                        if (block_i == block_j) continue;

                        assert(block_j < SizeType(block_pattern.size()));

                        if (block_pattern[block_j] == 0) {
                            assert((block_i + 1) < SizeType(row_ptr.size()));
                            row_ptr[block_i + 1]++;
                            block_pattern[block_j] = 1.0;
                        }
                    }
                }

                for (SizeType i = offset_i; i < offset_i + BlockSize; ++i) {
                    const SizeType row_begin = ia[i];
                    const SizeType row_end = ia[i + 1];

                    for (SizeType k = row_begin; k < row_end; ++k) {
                        SizeType j = ja[k];
                        SizeType block_j = j / BlockSize;

                        if (block_i == block_j) continue;

                        block_pattern[block_j] = 0;
                    }
                }
            }

            for (SizeType i = 0; i < n_blocks; ++i) {
                row_ptr[i + 1] += row_ptr[i];
            }

            colidx.resize(row_ptr[n_blocks], 0);
            std::vector<SizeType> count(n_blocks, 0);

            for (SizeType block_i = 0; block_i < n_blocks; block_i++) {
                const SizeType offset_i = block_i * BlockSize;

                for (SizeType i = offset_i; i < offset_i + BlockSize; ++i) {
                    const SizeType row_begin = ia[i];
                    const SizeType row_end = ia[i + 1];

                    for (SizeType k = row_begin; k < row_end; ++k) {
                        SizeType j = ja[k];
                        SizeType block_j = j / BlockSize;

                        if (block_i == block_j) {
                            //?
                        } else if (block_pattern[block_j] == 0) {
                            block_pattern[block_j] = 1.0;
                            SizeType idx_offset = row_ptr[block_i] + count[block_i];
                            colidx[idx_offset] = block_j;
                            ++count[block_i];
                        }
                    }
                }

                for (SizeType i = offset_i; i < offset_i + BlockSize; ++i) {
                    const SizeType row_begin = ia[i];
                    const SizeType row_end = ia[i + 1];

                    for (SizeType k = row_begin; k < row_end; ++k) {
                        SizeType j = ja[k];
                        SizeType block_j = j / BlockSize;
                        if (block_i != block_j) {
                            block_pattern[block_j] = 0;
                        }
                    }
                }

                if (sort_columns) {
                    std::sort(&colidx[row_ptr[block_i]], &colidx[row_ptr[block_i + 1] - 1]);
                }
            }

            values.resize(row_ptr[n_blocks] * BlockSize_2);

            update_split_diag(n, ia, ja, array, out, diag);
        }

        template <typename IA, typename JA, typename Array, typename Diag>
        static void update_split_diag(const SizeType n,
                                      const IA &ia,
                                      const JA &ja,
                                      const Array &array,
                                      CRSMatrix<OutScalarView, OutIndexView, BlockSize> &out,
                                      Diag &diag) {
            static const int BlockSize_2 = BlockSize * BlockSize;

            UTOPIA_UNUSED(n);

            const SizeType n_blocks = out.rows();
            assert(n == (BlockSize * n_blocks));

            auto &row_ptr = out.row_ptr();
            auto &colidx = out.colidx();
            auto &values = out.values();
            device::fill(0.0, values);

            for (SizeType block_i = 0; block_i < n_blocks; block_i++) {
                const SizeType offset_i = block_i * BlockSize;
                const SizeType block_row_begin = row_ptr[block_i];
                const SizeType block_row_end = row_ptr[block_i + 1];

                for (SizeType i = offset_i; i < offset_i + BlockSize; ++i) {
                    const SizeType row_begin = ia[i];
                    const SizeType row_end = ia[i + 1];
                    const SizeType sub_i = i - offset_i;

                    SizeType current_block = block_row_begin;

                    for (SizeType k = row_begin; k < row_end; ++k) {
                        const SizeType j = ja[k];
                        const Scalar value = array[k];

                        const SizeType block_j = j / BlockSize;
                        const SizeType sub_j = j - block_j * BlockSize;

                        if (block_i == block_j) {
                            diag[block_i * BlockSize_2 + sub_i * BlockSize + sub_j] = value;
                        } else {
                            SizeType k_block = current_block;
                            for (; k_block < block_row_end; ++k_block) {
                                if (colidx[k_block] == block_j) {
                                    current_block = k_block;
                                    break;
                                }
                            }

                            assert(colidx[k_block] == block_j);

                            values[k_block * BlockSize_2 + (sub_i * BlockSize) + sub_j] = value;
                        }
                    }
                }
            }
        }
    };

    template <class ScalarView, class IndexView, int BlockSize, class OutScalarView, class OutIndexView>
    void convert(const CRSMatrix<ScalarView, IndexView, 1> &mat,
                 CRSMatrix<OutScalarView, OutIndexView, BlockSize> &block_mat) {
        CRSToBlockCRS<ScalarView, IndexView, BlockSize, OutScalarView, OutIndexView>::apply(mat, block_mat);
    }

    template <class ScalarView, class IndexView, int BlockSize, class OutScalarView, class OutIndexView, class DiagView>
    void convert_split_diag(const CRSMatrix<ScalarView, IndexView, 1> &mat,
                            CRSMatrix<OutScalarView, OutIndexView, BlockSize> &block_mat,
                            DiagView &diag) {
        CRSToBlockCRS<ScalarView, IndexView, BlockSize, OutScalarView, OutIndexView>::apply_split_diag(
            mat, block_mat, diag);
    }

}  // namespace utopia

#endif  // UTOPIA_CRS_TO_BLOCK_CRS_HPP
