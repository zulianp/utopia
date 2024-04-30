#ifndef UTOPIA_CRS_TO_BLOCK_CRS_HPP
#define UTOPIA_CRS_TO_BLOCK_CRS_HPP

#include "utopia_Base.hpp"

#include "utopia_Algorithms.hpp"
#include "utopia_Traits.hpp"

#include "utopia_CRSMatrix.hpp"

#include <vector>

namespace utopia {

    template <typename SizeType, int BlockSize>
    class BlockCRSDefaultBlockAssignment {
    public:
        inline constexpr SizeType operator[](const SizeType i) const { return i / BlockSize; }
        inline constexpr int sub(const SizeType i) const { return i - (i / BlockSize) * BlockSize; }
    };

    template <class ScalarView,
              class IndexView,
              int BlockSize,
              class OutScalarView = ScalarView,
              class OutIndexView = IndexView>
    class CRSToBlockCRS {
    public:
        using Scalar = typename std::remove_const<typename Traits<ScalarView>::ValueType>::type;
        using SizeType = typename std::remove_const<typename Traits<IndexView>::ValueType>::type;

        using DefaultBlockAssignment = utopia::BlockCRSDefaultBlockAssignment<SizeType, BlockSize>;

        static void apply(const CRSMatrix<ScalarView, IndexView, 1> &in,
                          CRSMatrix<OutScalarView, OutIndexView, BlockSize> &out) {
            apply(in.rows(), in.row_ptr(), in.colidx(), in.values(), out);
        }

        template <typename RowBlockAssignemnt, typename ColBlockAssignemnt>
        static void apply(const CRSMatrix<ScalarView, IndexView, 1> &in,
                          const RowBlockAssignemnt &row_block_assignement,
                          const ColBlockAssignemnt &col_block_assignement,
                          CRSMatrix<OutScalarView, OutIndexView, BlockSize> &out) {
            apply(in.rows(), in.row_ptr(), in.colidx(), in.values(), row_block_assignement, col_block_assignement, out);
        }

        template <typename IA, typename JA, typename Array>
        static void apply(const SizeType n,
                          const IA &ia,
                          const JA &ja,
                          const Array &array,
                          CRSMatrix<OutScalarView, OutIndexView, BlockSize> &out) {
            apply(n, ia, ja, array, DefaultBlockAssignment(), DefaultBlockAssignment(), out);
        }

        template <typename IA, typename JA, typename Array, typename RowBlockAssignemnt, typename ColBlockAssignemnt>
        static void apply(const SizeType n,
                          const IA &ia,
                          const JA &ja,
                          const Array &array,
                          const RowBlockAssignemnt &row_block_assignement,
                          const ColBlockAssignemnt &col_block_assignement,
                          CRSMatrix<OutScalarView, OutIndexView, BlockSize> &out) {
            static const int BlockSize_2 = BlockSize * BlockSize;
            auto n_blocks = row_block_assignement[n - 1] + 1;

            // assert(n_blocks * BlockSize == n);
            assert(n > 0);

            auto n_nnz = ia[n];
            SizeType max_j = 0;

            for (SizeType i = 0; i < n_nnz; ++i) {
                max_j = std::max(max_j, ja[i]);
            }

            max_j = col_block_assignement[max_j];

            std::vector<SizeType> block_lut(max_j + 1, -1);
            std::vector<SizeType> current_row(max_j + 1, 0);

            // FIXME
            SizeType cols = std::max(n_blocks, max_j + 1);
            out.set_cols(cols);

            auto &row_ptr = out.row_ptr();
            auto &colidx = out.colidx();
            auto &values = out.values();

            row_ptr.resize(n_blocks + 1);
            device::fill(SizeType(0), row_ptr);

            SizeType prev_block_i = 0;

            // Detect number of blocks per block row
            for (SizeType i = 0; i < n; ++i) {
                const SizeType row_begin = ia[i];
                const SizeType row_end = ia[i + 1];
                // const SizeType block_i = i / BlockSize;
                const SizeType block_i = row_block_assignement[i];

                if (prev_block_i != block_i) {
                    // Clean-up boolean flags
                    const SizeType n_blocks_in_row = row_ptr[prev_block_i + 1];
                    for (SizeType k = 0; k < n_blocks_in_row; ++k) {
                        assert(k < SizeType(current_row.size()));
                        SizeType block_j = current_row[k];
                        assert(block_j >= 0);
                        assert(block_j < (max_j + 1));
                        block_lut[block_j] = -1;
                    }

                    prev_block_i = block_i;
                }

                for (SizeType k = row_begin; k < row_end; ++k) {
                    const SizeType j = ja[k];
                    // const SizeType block_j = j / BlockSize;
                    const SizeType block_j = col_block_assignement[j];

                    assert(block_j >= 0);
                    assert(block_j < SizeType(block_lut.size()));

                    if (block_lut[block_j] == -1) {
                        assert((block_i + 1) < SizeType(row_ptr.size()));

                        // Store current column index
                        current_row[row_ptr[block_i + 1]] = block_j;

                        // Flag current column as encountered
                        block_lut[block_j] = 1;  // use as boolean

                        // Increase number of columns in block row
                        row_ptr[block_i + 1]++;

                        // Check that we do not overflow
                        assert(row_ptr[block_i + 1] <= cols);
                    }
                }
            }

            // Compute row pointer
            for (SizeType i = 0; i < n_blocks; ++i) {
                row_ptr[i + 1] += row_ptr[i];
                assert((row_ptr[i + 1] - row_ptr[i]) <= cols);
                assert(row_ptr[i + 1] >= 0);
            }

#ifndef NDEBUG
            for (SizeType i = 0; i < n_blocks; ++i) {
                assert(row_ptr[i + 1] >= row_ptr[i]);
                assert((row_ptr[i + 1] - row_ptr[i]) <= cols);
            }
#endif  // NDEBUG

            // Allocate buffers
            colidx.resize(row_ptr[n_blocks], 0);

            // FIXME
            device::fill(-1, block_lut);
            values.resize(row_ptr[n_blocks] * BlockSize_2);
            device::fill(0., values);

            prev_block_i = 0;
            SizeType count = 0;
            for (SizeType i = 0; i < n; ++i) {
                const SizeType row_begin = ia[i];
                const SizeType row_end = ia[i + 1];
                // const SizeType block_i = i / BlockSize;
                const SizeType block_i = row_block_assignement[i];
                const SizeType sub_i = row_block_assignement.sub(i);

                assert(sub_i < BlockSize);
                assert(sub_i >= 0);

                if (prev_block_i != block_i) {
                    // Clean-up boolean flags
                    const SizeType n_blocks_in_row = row_ptr[prev_block_i + 1] - row_ptr[prev_block_i];
                    assert(n_blocks_in_row <= cols);
                    for (SizeType k = 0; k < n_blocks_in_row; ++k) {
                        assert(k < SizeType(current_row.size()));
                        SizeType block_j = current_row[k];
                        assert(block_j >= 0);
                        assert(block_j < SizeType(block_lut.size()));
                        block_lut[block_j] = -1;
                    }

                    prev_block_i = block_i;
                    count = 0;
                }

                for (SizeType k = row_begin; k < row_end; ++k) {
                    const SizeType j = ja[k];
                    const Scalar value = array[k];
                    // const SizeType block_j = j / BlockSize;
                    const SizeType block_j = col_block_assignement[j];
                    const SizeType sub_j = col_block_assignement.sub(j);

                    assert(sub_j < BlockSize);
                    assert(sub_j >= 0);

                    assert(block_j < SizeType(block_lut.size()));

                    if (block_lut[block_j] == -1) {
                        assert((block_i + 1) < SizeType(row_ptr.size()));

                        // Store current column index
                        assert(count < SizeType(current_row.size()));
                        current_row[count] = block_j;

                        // Set column index in new slot
                        SizeType idx_offset = row_ptr[block_i] + count;

                        // Store current column index
                        block_lut[block_j] = idx_offset;

                        colidx[idx_offset] = block_j;

                        // Increase counter
                        ++count;
                    }

                    // Set values in block
                    assert(block_lut[block_j] != -1);
                    SizeType val_idx = block_lut[block_j] * BlockSize_2 + (sub_i * BlockSize) + sub_j;
                    assert(val_idx < SizeType(values.size()));
                    values[val_idx] = value;
                }
            }

            // if (sort_columns) {
            //     for (SizeType block_i = 0; block_i < n_blocks; block_i++) {
            //         if (row_ptr[block_i] < row_ptr[block_i + 1]) {
            //             std::sort(&colidx[0] + row_ptr[block_i], &colidx[0] + row_ptr[block_i + 1]);
            //         }
            //     }
            // }

            // values.resize(row_ptr[n_blocks] * BlockSize_2);

            // update(n, ia, ja, array, row_block_assignement, col_block_assignement, out);

            assert(out.is_valid());
        }

        static void update(const CRSMatrix<ScalarView, IndexView, 1> &in,
                           CRSMatrix<OutScalarView, OutIndexView, BlockSize> &out) {
            update(in.rows(),
                   in.row_ptr(),
                   in.colidx(),
                   in.values(),
                   DefaultBlockAssignment(),
                   DefaultBlockAssignment(),
                   out);
        }

        // template <typename IA, typename JA, typename Array>
        // static void update(const SizeType n,
        //                    const IA &ia,
        //                    const JA &ja,
        //                    const Array &array,
        //                    CRSMatrix<OutScalarView, OutIndexView, BlockSize> &out) {

        // }

        template <typename IA, typename JA, typename Array, class RowBlockAssignemnt, class ColBlockAssignemnt>
        static void update(const SizeType n,
                           const IA &ia,
                           const JA &ja,
                           const Array &array,
                           const RowBlockAssignemnt &,
                           const ColBlockAssignemnt &,
                           CRSMatrix<OutScalarView, OutIndexView, BlockSize> &out) {
            static const int BlockSize_2 = BlockSize * BlockSize;

            UTOPIA_UNUSED(n);

            const SizeType n_blocks = out.rows();
            assert(n == (BlockSize * n_blocks));

            const auto &row_ptr = out.row_ptr();
            const auto &colidx = out.colidx();
            auto &values = out.values();
            device::fill(0.0, values);

            // for (SizeType i = 0; i < n; ++i) {
            //     const SizeType row_begin = ia[i];
            //     const SizeType row_end = ia[i + 1];
            //     // const SizeType block_i = i / BlockSize;
            //     const SizeType block_i = row_block_assignement[i];
            //     const SizeType sub_i = row_block_assignement.sub(i);

            //     for (SizeType k = row_begin; k < row_end; ++k) {
            //         const SizeType j = ja[k];
            //         // const SizeType block_j = j / BlockSize;
            //         const SizeType block_j = col_block_assignement[j];
            //         const SizeType sub_j = col_block_assignement.sub(j);

            //         values[k_block * BlockSize_2 + (sub_i * BlockSize) + sub_j] = value;
            //     }
            // }

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

            const auto &row_ptr = out.row_ptr();
            const auto &colidx = out.colidx();
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

        template <typename Diag>
        static void update_split_diag(const CRSMatrix<ScalarView, IndexView, 1> &in,
                                      CRSMatrix<OutScalarView, OutIndexView, BlockSize> &out,
                                      Diag &diag) {
            update_split_diag(in.rows(), in.row_ptr(), in.colidx(), in.values(), out, diag);
        }
    };

    template <class ScalarView, class IndexView, int BlockSize, class OutScalarView, class OutIndexView>
    void convert(const CRSMatrix<ScalarView, IndexView, 1> &mat,
                 CRSMatrix<OutScalarView, OutIndexView, BlockSize> &block_mat) {
        CRSToBlockCRS<ScalarView, IndexView, BlockSize, OutScalarView, OutIndexView>::apply(mat, block_mat);
    }

    template <class ScalarView,
              class IndexView,
              int BlockSize,
              class RowBlockAssignement,
              class ColBlockAssignement,
              class OutScalarView,
              class OutIndexView>
    void convert(const CRSMatrix<ScalarView, IndexView, 1> &mat,
                 const RowBlockAssignement &row_block_assignement,
                 const ColBlockAssignement &col_block_assignement,
                 CRSMatrix<OutScalarView, OutIndexView, BlockSize> &block_mat) {
        CRSToBlockCRS<ScalarView, IndexView, BlockSize, OutScalarView, OutIndexView>::apply(
            mat, row_block_assignement, col_block_assignement, block_mat);
    }

    template <class ScalarView, class IndexView, int BlockSize, class OutScalarView, class OutIndexView, class DiagView>
    void convert_split_diag(const CRSMatrix<ScalarView, IndexView, 1> &mat,
                            CRSMatrix<OutScalarView, OutIndexView, BlockSize> &block_mat,
                            DiagView &diag) {
        CRSToBlockCRS<ScalarView, IndexView, BlockSize, OutScalarView, OutIndexView>::apply_split_diag(
            mat, block_mat, diag);
    }

    template <class ScalarView, class IndexView, int BlockSize, class OutScalarView, class OutIndexView, class DiagView>
    void convert_split_diag_update(const CRSMatrix<ScalarView, IndexView, 1> &mat,
                                   CRSMatrix<OutScalarView, OutIndexView, BlockSize> &block_mat,
                                   DiagView &diag) {
        CRSToBlockCRS<ScalarView, IndexView, BlockSize, OutScalarView, OutIndexView>::update_split_diag(
            mat, block_mat, diag);
    }

}  // namespace utopia

#endif  // UTOPIA_CRS_TO_BLOCK_CRS_HPP
