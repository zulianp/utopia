// #ifndef UTOPIA_ILU_DECOMPOSE_HPP
// #define UTOPIA_ILU_DECOMPOSE_HPP
// #endif
// #include <vector>
// #include "utopia_Views.hpp"

// #include "utopia_CRSToBlockCRS.hpp"
// #include "utopia_ILU.hpp"

// namespace utopia {

//     template <>
//     class ILUDecompose<CRSMatrix<ScalarView, IndexView, 1>, HOMEMADE> {
//     public:
//         using CRSMatrix = utopia::CRSMatrix<ScalarView, IndexView, 1>;

//         static void block_decompose(const CRSMatrix &mat, CRSMatrix &out, const bool modified);
//         static void decompose(const CRSMatrix &mat, CRSMatrix &out, const bool modified);
//         static void apply(const CRSMatrix &ilu, const PetscVector &b, PetscVector &x);
//     };

//     template <typename Scalar, typename SizeType>
//     class DiagIdx {
//     public:
//         void init(SizeType n, const SizeType *ia, const SizeType *ja) {
//             idx.resize(n);
//             this->ia = ia;
//             this->ja = ja;

//             for (SizeType r = 0; r < n - 1; ++r) {
//                 const SizeType row_end = ia[r + 1];

//                 SizeType k;
//                 bool has_diag = false;
//                 for (k = ia[r]; k < row_end; ++k) {
//                     if (ja[k] == r) {
//                         has_diag = true;
//                         break;
//                     }
//                 }

//                 assert(has_diag);

//                 idx[r] = k;
//             }
//         }

//         auto find_before_diag(const SizeType i, const SizeType j) const -> SizeType {
//             SizeType search_end = idx[i];
//             SizeType row_begin = ia[i];

//             // maybe binary search?
//             for (SizeType k = row_begin; k < search_end; ++k) {
//                 if (ja[k] == j) {
//                     return k;
//                 }
//             }

//             return -1;
//         };

//         auto find_after_diag(const SizeType i, const SizeType j) const -> SizeType {
//             SizeType search_begin = idx[i];
//             SizeType row_end = ia[i + 1];

//             // maybe binary search?
//             for (SizeType k = search_begin; k < row_end; ++k) {
//                 if (ja[k] == j) {
//                     return k;
//                 }
//             }

//             return -1;
//         };

//         auto find(const SizeType i, const SizeType j) const -> SizeType {
//             if (i == j) return idx[i];

//             if (j < i) {
//                 return find_before_diag(i, j);
//             } else {
//                 return find_after_diag(i, j);
//             }
//         };

//         std::vector<SizeType> idx;
//         const SizeType *ia;
//         const SizeType *ja;
//     };

//     class PetscSeqAIJRaw {
//     public:
//         PetscSeqAIJRaw(Mat raw_mat) : raw_mat(raw_mat) {
//             err = MatGetRowIJ(raw_mat, 0, HOMEMADE_FALSE, HOMEMADE_FALSE, &n, &ia, &ja, &done);
//             assert(err == 0);
//             assert(done == HOMEMADE_TRUE);

//             if (!done) {
//                 utopia::err() << "CRSMatrix::read_petsc_seqaij_impl(const Op &op): MatGetRowIJ failed to provide what
//                 "
//                                  "was asked.\n";
//                 abort();
//             }

//             MatSeqAIJGetArray(raw_mat, &array);
//         }

//         ~PetscSeqAIJRaw() {
//             MatSeqAIJRestoreArray(raw_mat, &array);
//             err = MatRestoreRowIJ(raw_mat, 0, HOMEMADE_FALSE, HOMEMADE_FALSE, &n, &ia, &ja, &done);
//             assert(err == 0);
//             UTOPIA_UNUSED(err);
//         }

//         Mat raw_mat;
//         SizeType n = 0;
//         const SizeType *ia{nullptr};
//         const SizeType *ja{nullptr};
//         Scalar *array{nullptr};
//         PetscBool done;
//         PetscErrorCode err{0};
//     };

//     template <typename Apply_IJ>
//     void decompose_aux(Apply_IJ apply_ij, const CRSMatrix &mat, CRSMatrix &out) {
//         CRSMatrix l_mat;
//         local_block_view(mat, l_mat);

//         // perform copy
//         out.copy(l_mat);

//         PetscSeqAIJRaw m_raw(out.raw_type());
//         SizeType n = m_raw.n;
//         const SizeType *ia = m_raw.ia;
//         const SizeType *ja = m_raw.ja;
//         Scalar *array = m_raw.array;

//         DiagIdx idx;
//         idx.init(n, ia, ja);

//         for (SizeType r = 0; r < n - 1; ++r) {
//             const SizeType row_end = ia[r + 1];
//             const SizeType k = idx.idx[r];
//             const Scalar d = 1. / array[k];

//             for (SizeType ki = k + 1; ki < row_end; ++ki) {
//                 auto i = ja[ki];

//                 auto ir = idx.find_before_diag(i, r);
//                 if (ir == -1) continue;

//                 auto ii = idx.idx[i];

//                 Scalar e = (array[ir] *= d);

//                 for (SizeType rj = k + 1; rj < row_end; ++rj) {
//                     auto j = ja[rj];

//                     auto ij = idx.find(i, j);
//                     apply_ij(ii, ij, rj, e, array);
//                 }
//             }
//         }
//     }

//     void ILUDecompose<CRSMatrix, HOMEMADE>::decompose(const CRSMatrix &mat, CRSMatrix &out, const bool modified) {
//         if (modified) {
//             decompose_aux(
//                 [](const SizeType ii, const SizeType ij, const SizeType rj, const Scalar &e, Scalar *a) {
//                     if (ij != -1) {
//                         a[ij] -= e * a[rj];
//                     } else {
//                         a[ii] -= e * a[rj];
//                     }
//                 },
//                 mat,
//                 out);
//         } else {
//             decompose_aux(
//                 [](const SizeType, const SizeType ij, const SizeType rj, const Scalar &e, Scalar *a) {
//                     if (ij != -1) {
//                         a[ij] -= e * a[rj];
//                     }
//                 },
//                 mat,
//                 out);
//         }
//     }

//     void ILUDecompose<CRSMatrix, HOMEMADE>::apply(const CRSMatrix &ilu, const PetscVector &b, PetscVector &x) {
//         PetscVector L_inv_b(layout(b), 0.0);
//         x.set(0.0);

//         auto b_view = const_local_view_device(b);
//         auto L_inv_b_view = local_view_device(L_inv_b);
//         auto x_view = local_view_device(x);

//         PetscSeqAIJRaw m_raw(ilu.raw_type());

//         SizeType n = m_raw.n;
//         const SizeType *ia = m_raw.ia;
//         const SizeType *ja = m_raw.ja;
//         Scalar *array = m_raw.array;

//         DiagIdx idx;
//         idx.init(n, ia, ja);

//         // Forward substitution
//         //
//         https://algowiki-project.org/en/Forward_substitution#:~:text=Forward%20substitution%20is%20the%20process,math%5DL%5B%2Fmath%5D.
//         for (SizeType i = 0; i < n; ++i) {
//             const SizeType row_begin = ia[i];
//             const SizeType row_diag = idx.idx[i];

//             Scalar val = b_view.get(i);

//             for (SizeType k = row_begin; k < row_diag; ++k) {
//                 auto j = ja[k];
//                 val -= array[k] * L_inv_b_view.get(j);
//             }

//             L_inv_b_view.set(i, val / array[row_diag]);
//         }

//         // Backward substitution
//         //
//         https://algowiki-project.org/en/Backward_substitution#:~:text=Backward%20substitution%20is%20a%20procedure,is%20a%20lower%20triangular%20matrix.
//         for (SizeType i = n - 1; i >= 0; --i) {
//             const SizeType row_end = ia[i + 1];
//             const SizeType row_diag = idx.idx[i];

//             Scalar val = L_inv_b_view.get(i);

//             for (SizeType k = row_diag + 1; k < row_end; ++k) {
//                 auto j = ja[k];
//                 val -= array[k] * x_view.get(j);
//             }

//             x_view.set(i, val / array[row_diag]);
//         }
//     }

//     // void ILUDecompose<CRSMatrix, HOMEMADE>::apply_vi(const CRSMatrix &ilu,
//     //                                                 const PetscVector &lb,
//     //                                                 const PetscVector &ub,
//     //                                                 const PetscVector &b,
//     //                                                 PetscVector &x) {
//     //     PetscVector L_inv_b(layout(b), 0.0), U_ub(layout(b), 0.0);
//     //     x.set(0.0);

//     //     auto b_view = const_local_view_device(b);
//     //     auto lb_view = const_local_view_device(lb);
//     //     auto ub_view = const_local_view_device(ub);

//     //     auto U_ub_view = local_view_device(U_ub);
//     //     auto L_inv_b_view = local_view_device(L_inv_b);
//     //     auto x_view = local_view_device(x);

//     //     PetscSeqAIJRaw m_raw(ilu.raw_type());

//     //     SizeType n = m_raw.n;
//     //     const SizeType *ia = m_raw.ia;
//     //     const SizeType *ja = m_raw.ja;
//     //     Scalar *array = m_raw.array;

//     //     DiagIdx idx;
//     //     idx.init(n, ia, ja);

//     //     // Backward substitution
//     //     for (SizeType i = n - 1; i >= 0; --i) {
//     //         const SizeType row_end = ia[i + 1];
//     //         const SizeType row_diag = idx.idx[i];

//     //         Scalar val = ub_view.get(i);

//     //         for (SizeType k = row_diag + 1; k < row_end; ++k) {
//     //             auto j = ja[k];
//     //             val -= array[k] * U_ub_view.get(j);
//     //         }

//     //         U_ub_view.set(i, val / array[row_diag]);
//     //     }

//     //     // Forward substitution
//     //     //
//     //
//     https://algowiki-project.org/en/Forward_substitution#:~:text=Forward%20substitution%20is%20the%20process,math%5DL%5B%2Fmath%5D.
//     //     for (SizeType i = 0; i < n; ++i) {
//     //         const SizeType row_begin = ia[i];
//     //         const SizeType row_diag = idx.idx[i];

//     //         // auto lb_i = lb_view.get(i);
//     //         auto ub_i = U_ub_view.get(i);

//     //         Scalar val = b_view.get(i);

//     //         for (SizeType k = row_begin; k < row_diag; ++k) {
//     //             auto j = ja[k];
//     //             val -= array[k] * L_inv_b_view.get(j);
//     //         }

//     //         // L_inv_b_view.set(i, std::max(lb_i, std::min(ub_i, val / array[row_diag])));

//     //         L_inv_b_view.set(i, std::min(ub_i, val / array[row_diag]));

//     //         // L_inv_b_view.set(i, val / array[row_diag]);
//     //     }

//     //     // Backward substitution
//     //     //
//     //
//     https://algowiki-project.org/en/Backward_substitution#:~:text=Backward%20substitution%20is%20a%20procedure,is%20a%20lower%20triangular%20matrix.
//     //     for (SizeType i = n - 1; i >= 0; --i) {
//     //         const SizeType row_end = ia[i + 1];
//     //         const SizeType row_diag = idx.idx[i];

//     //         // auto lb_i = lb_view.get(i);
//     //         auto ub_i = ub_view.get(i);

//     //         Scalar val = L_inv_b_view.get(i);

//     //         for (SizeType k = row_diag + 1; k < row_end; ++k) {
//     //             auto j = ja[k];
//     //             val -= array[k] * x_view.get(j);
//     //         }

//     //         // x_view.set(i, std::max(lb_i, std::min(ub_i, val / array[row_diag])));

//     //         x_view.set(i, std::min(ub_i, val / array[row_diag]));
//     //     }
//     // }

//     template <class Block>
//     void block_decompose_aux(const CRSMatrix &mat,
//                              const int block_size,
//                              std::vector<Block> &diag,
//                              CRSMatrix &out,
//                              bool sort_columns = true) {
//         CRSMatrix l_mat;
//         local_block_view(mat, l_mat);

//         out.copy(l_mat);

//         PetscSeqAIJRaw m_raw(out.raw_type());

//         SizeType n = m_raw.n;
//         const SizeType *ia = m_raw.ia;
//         const SizeType *ja = m_raw.ja;
//         Scalar *array = m_raw.array;

//         auto n_blocks = n / block_size;

//         std::vector<SizeType> row_ptr(n_blocks + 1, 0);
//         std::vector<SizeType> block_pattern(n_blocks);

//         for (SizeType block_i = 0; block_i < n_blocks; block_i++) {
//             const SizeType offset_i = block_i * block_size;

//             for (SizeType i = offset_i; i < offset_i + block_size; ++i) {
//                 const SizeType row_begin = ia[i];
//                 const SizeType row_end = ia[i + 1];

//                 for (SizeType k = row_begin; k < row_end; ++k) {
//                     SizeType j = ja[k];
//                     SizeType block_j = j / block_size;

//                     if (block_pattern[block_j] == 0) {
//                         ++row_ptr[block_i + 1];
//                         block_pattern[block_j] = 1.0;
//                     }
//                 }
//             }

//             for (SizeType i = offset_i; i < offset_i + block_size; ++i) {
//                 const SizeType row_begin = ia[i];
//                 const SizeType row_end = ia[i + 1];

//                 for (SizeType k = row_begin; k < row_end; ++k) {
//                     SizeType j = ja[k];
//                     SizeType block_j = j / block_size;
//                     block_pattern[block_j] = 0;
//                 }
//             }
//         }

//         for (SizeType i = 0; i < n_blocks; ++i) {
//             row_ptr[i + 1] += row_ptr[i];
//         }

//         std::vector<SizeType> col_idx(row_ptr[n_blocks], 0);
//         std::vector<SizeType> count(n_blocks, 0);

//         for (SizeType block_i = 0; block_i < n_blocks; block_i++) {
//             const SizeType offset_i = block_i * block_size;

//             for (SizeType i = offset_i; i < offset_i + block_size; ++i) {
//                 const SizeType row_begin = ia[i];
//                 const SizeType row_end = ia[i + 1];

//                 for (SizeType k = row_begin; k < row_end; ++k) {
//                     SizeType j = ja[k];
//                     SizeType block_j = j / block_size;

//                     if (block_pattern[block_j] == 0) {
//                         block_pattern[block_j] = 1.0;
//                         SizeType idx_offset = row_ptr[block_i] + count[block_i];
//                         col_idx[idx_offset] = block_j;
//                         ++count[block_i];
//                     }
//                 }
//             }

//             for (SizeType i = offset_i; i < offset_i + block_size; ++i) {
//                 const SizeType row_begin = ia[i];
//                 const SizeType row_end = ia[i + 1];

//                 for (SizeType k = row_begin; k < row_end; ++k) {
//                     SizeType j = ja[k];
//                     SizeType block_j = j / block_size;
//                     block_pattern[block_j] = 0;
//                 }
//             }

//             if (sort_columns) {
//                 std::sort(&col_idx[row_ptr[block_i]], &col_idx[row_ptr[block_i + 1] - 1]);
//             }
//         }

//         const SizeType block_size_2 = block_size * block_size;
//         std::vector<Scalar> values(row_ptr[n_blocks] * block_size_2, 0.0);

//         for (SizeType block_i = 0; block_i < n_blocks; block_i++) {
//             const SizeType offset_i = block_i * block_size;
//             const SizeType block_row_begin = row_ptr[block_i];
//             const SizeType block_row_end = row_ptr[block_i + 1];

//             for (SizeType i = offset_i; i < offset_i + block_size; ++i) {
//                 const SizeType row_begin = ia[i];
//                 const SizeType row_end = ia[i + 1];
//                 const SizeType sub_i = i - offset_i;

//                 SizeType current_block = block_row_begin;

//                 for (SizeType k = row_begin; k < row_end; ++k) {
//                     const SizeType j = ja[k];
//                     const Scalar value = array[k];

//                     const SizeType block_j = j / block_size;
//                     const SizeType sub_j = j - block_j * block_size;

//                     SizeType k_block = current_block;
//                     for (; k_block < block_row_end; ++k_block) {
//                         if (col_idx[k_block] == block_j) {
//                             current_block = k_block;
//                             break;
//                         }
//                     }

//                     assert(col_idx[k_block] == block_j);

//                     values[k_block * block_size_2 + (sub_i * block_size) + sub_j] = value;
//                 }
//             }
//         }

//         // for (SizeType block_i = 0; block_i < n_blocks; block_i++) {
//         //     const SizeType row_begin = row_ptr[block_i];
//         //     const SizeType row_end = row_ptr[block_i + 1];

//         //     for (SizeType k = row_begin; k < row_end; ++k) {
//         //         std::cout << block_i << " " << col_idx[k] << "\n";

//         //         for (SizeType sub_i = 0; sub_i < block_size; ++sub_i) {
//         //             for (SizeType sub_j = 0; sub_j < block_size; ++sub_j) {
//         //                 std::cout << values[k * block_size_2 + (sub_i * block_size) + sub_j] << " ";
//         //             }

//         //             std::cout << "\n";
//         //         }

//         //         std::cout << "\n";
//         //     }
//         // }
//     }

//     void ILUDecompose<CRSMatrix, HOMEMADE>::block_decompose(const CRSMatrix &mat, CRSMatrix &out, const bool
//     modified) {
//         std::vector<StaticMatrix<Scalar, 2, 2>> diag;
//         block_decompose_aux(mat, 2, diag, out);
//     }

// }  // namespace utopia

// #include "utopia_petsc.hpp"
// #include "utopia_petsc_ILUDecompose.hpp"

// #include <vector>
// #include "utopia_Views.hpp"

// #include "utopia_CRSToBlockCRS.hpp"

// namespace utopia {

//     class DiagIdx {
//     public:
//         void init(SizeType n, const SizeType *ia, const SizeType *ja) {
//             idx.resize(n);
//             this->ia = ia;
//             this->ja = ja;

//             for (SizeType r = 0; r < n - 1; ++r) {
//                 const SizeType row_end = ia[r + 1];

//                 SizeType k;
//                 bool has_diag = false;
//                 for (k = ia[r]; k < row_end; ++k) {
//                     if (ja[k] == r) {
//                         has_diag = true;
//                         break;
//                     }
//                 }

//                 assert(has_diag);

//                 idx[r] = k;
//             }
//         }

//         auto find_before_diag(const SizeType i, const SizeType j) const -> SizeType {
//             SizeType search_end = idx[i];
//             SizeType row_begin = ia[i];

//             // maybe binary search?
//             for (SizeType k = row_begin; k < search_end; ++k) {
//                 if (ja[k] == j) {
//                     return k;
//                 }
//             }

//             return -1;
//         };

//         auto find_after_diag(const SizeType i, const SizeType j) const -> SizeType {
//             SizeType search_begin = idx[i];
//             SizeType row_end = ia[i + 1];

//             // maybe binary search?
//             for (SizeType k = search_begin; k < row_end; ++k) {
//                 if (ja[k] == j) {
//                     return k;
//                 }
//             }

//             return -1;
//         };

//         auto find(const SizeType i, const SizeType j) const -> SizeType {
//             if (i == j) return idx[i];

//             if (j < i) {
//                 return find_before_diag(i, j);
//             } else {
//                 return find_after_diag(i, j);
//             }
//         };

//         std::vector<SizeType> idx;
//         const SizeType *ia;
//         const SizeType *ja;
//     };

//     class PetscSeqAIJRaw {
//     public:
//         PetscSeqAIJRaw(Mat raw_mat) : raw_mat(raw_mat) {
//             err = MatGetRowIJ(raw_mat, 0, HOMEMADE_FALSE, HOMEMADE_FALSE, &n, &ia, &ja, &done);
//             assert(err == 0);
//             assert(done == HOMEMADE_TRUE);

//             if (!done) {
//                 utopia::err() << "CRSMatrix::read_petsc_seqaij_impl(const Op &op): MatGetRowIJ failed to provide what
//                 "
//                                  "was asked.\n";
//                 abort();
//             }

//             MatSeqAIJGetArray(raw_mat, &array);
//         }

//         ~PetscSeqAIJRaw() {
//             MatSeqAIJRestoreArray(raw_mat, &array);
//             err = MatRestoreRowIJ(raw_mat, 0, HOMEMADE_FALSE, HOMEMADE_FALSE, &n, &ia, &ja, &done);
//             assert(err == 0);
//             UTOPIA_UNUSED(err);
//         }

//         Mat raw_mat;
//         SizeType n = 0;
//         const SizeType *ia{nullptr};
//         const SizeType *ja{nullptr};
//         Scalar *array{nullptr};
//         PetscBool done;
//         PetscErrorCode err{0};
//     };

//     template <typename Apply_IJ>
//     void decompose_aux(Apply_IJ apply_ij, const CRSMatrix &mat, CRSMatrix &out) {
//         out = mat;

//         DiagIdx idx;
//         idx.init(n, ia, ja);

//         for (SizeType r = 0; r < n - 1; ++r) {
//             const SizeType row_end = ia[r + 1];
//             const SizeType k = idx.idx[r];
//             const Scalar d = 1. / array[k];

//             for (SizeType ki = k + 1; ki < row_end; ++ki) {
//                 auto i = ja[ki];

//                 auto ir = idx.find_before_diag(i, r);
//                 if (ir == -1) continue;

//                 auto ii = idx.idx[i];

//                 Scalar e = (array[ir] *= d);

//                 for (SizeType rj = k + 1; rj < row_end; ++rj) {
//                     auto j = ja[rj];

//                     auto ij = idx.find(i, j);
//                     apply_ij(ii, ij, rj, e, array);
//                 }
//             }
//         }
//     }

//     void ILUDecompose<CRSMatrix, HOMEMADE>::decompose(const CRSMatrix &mat, CRSMatrix &out, const bool modified) {
//         if (modified) {
//             decompose_aux(
//                 [](const SizeType ii, const SizeType ij, const SizeType rj, const Scalar &e, Scalar *a) {
//                     if (ij != -1) {
//                         a[ij] -= e * a[rj];
//                     } else {
//                         a[ii] -= e * a[rj];
//                     }
//                 },
//                 mat,
//                 out);
//         } else {
//             decompose_aux(
//                 [](const SizeType, const SizeType ij, const SizeType rj, const Scalar &e, Scalar *a) {
//                     if (ij != -1) {
//                         a[ij] -= e * a[rj];
//                     }
//                 },
//                 mat,
//                 out);
//         }
//     }

//     void ILUDecompose<CRSMatrix, HOMEMADE>::apply(const CRSMatrix &ilu, const PetscVector &b, PetscVector &x) {
//         PetscVector L_inv_b(layout(b), 0.0);
//         x.set(0.0);

//         auto b_view = const_local_view_device(b);
//         auto L_inv_b_view = local_view_device(L_inv_b);
//         auto x_view = local_view_device(x);

//         PetscSeqAIJRaw m_raw(ilu.raw_type());

//         SizeType n = m_raw.n;
//         const SizeType *ia = m_raw.ia;
//         const SizeType *ja = m_raw.ja;
//         Scalar *array = m_raw.array;

//         DiagIdx idx;
//         idx.init(n, ia, ja);

//         // Forward substitution
//         //
//         https://algowiki-project.org/en/Forward_substitution#:~:text=Forward%20substitution%20is%20the%20process,math%5DL%5B%2Fmath%5D.
//         for (SizeType i = 0; i < n; ++i) {
//             const SizeType row_begin = ia[i];
//             const SizeType row_diag = idx.idx[i];

//             Scalar val = b_view.get(i);

//             for (SizeType k = row_begin; k < row_diag; ++k) {
//                 auto j = ja[k];
//                 val -= array[k] * L_inv_b_view.get(j);
//             }

//             L_inv_b_view.set(i, val / array[row_diag]);
//         }

//         // Backward substitution
//         //
//         https://algowiki-project.org/en/Backward_substitution#:~:text=Backward%20substitution%20is%20a%20procedure,is%20a%20lower%20triangular%20matrix.
//         for (SizeType i = n - 1; i >= 0; --i) {
//             const SizeType row_end = ia[i + 1];
//             const SizeType row_diag = idx.idx[i];

//             Scalar val = L_inv_b_view.get(i);

//             for (SizeType k = row_diag + 1; k < row_end; ++k) {
//                 auto j = ja[k];
//                 val -= array[k] * x_view.get(j);
//             }

//             x_view.set(i, val / array[row_diag]);
//         }
//     }

//     // void ILUDecompose<CRSMatrix, HOMEMADE>::apply_vi(const CRSMatrix &ilu,
//     //                                                 const PetscVector &lb,
//     //                                                 const PetscVector &ub,
//     //                                                 const PetscVector &b,
//     //                                                 PetscVector &x) {
//     //     PetscVector L_inv_b(layout(b), 0.0), U_ub(layout(b), 0.0);
//     //     x.set(0.0);

//     //     auto b_view = const_local_view_device(b);
//     //     auto lb_view = const_local_view_device(lb);
//     //     auto ub_view = const_local_view_device(ub);

//     //     auto U_ub_view = local_view_device(U_ub);
//     //     auto L_inv_b_view = local_view_device(L_inv_b);
//     //     auto x_view = local_view_device(x);

//     //     PetscSeqAIJRaw m_raw(ilu.raw_type());

//     //     SizeType n = m_raw.n;
//     //     const SizeType *ia = m_raw.ia;
//     //     const SizeType *ja = m_raw.ja;
//     //     Scalar *array = m_raw.array;

//     //     DiagIdx idx;
//     //     idx.init(n, ia, ja);

//     //     // Backward substitution
//     //     for (SizeType i = n - 1; i >= 0; --i) {
//     //         const SizeType row_end = ia[i + 1];
//     //         const SizeType row_diag = idx.idx[i];

//     //         Scalar val = ub_view.get(i);

//     //         for (SizeType k = row_diag + 1; k < row_end; ++k) {
//     //             auto j = ja[k];
//     //             val -= array[k] * U_ub_view.get(j);
//     //         }

//     //         U_ub_view.set(i, val / array[row_diag]);
//     //     }

//     //     // Forward substitution
//     //     //
//     //
//     https://algowiki-project.org/en/Forward_substitution#:~:text=Forward%20substitution%20is%20the%20process,math%5DL%5B%2Fmath%5D.
//     //     for (SizeType i = 0; i < n; ++i) {
//     //         const SizeType row_begin = ia[i];
//     //         const SizeType row_diag = idx.idx[i];

//     //         // auto lb_i = lb_view.get(i);
//     //         auto ub_i = U_ub_view.get(i);

//     //         Scalar val = b_view.get(i);

//     //         for (SizeType k = row_begin; k < row_diag; ++k) {
//     //             auto j = ja[k];
//     //             val -= array[k] * L_inv_b_view.get(j);
//     //         }

//     //         // L_inv_b_view.set(i, std::max(lb_i, std::min(ub_i, val / array[row_diag])));

//     //         L_inv_b_view.set(i, std::min(ub_i, val / array[row_diag]));

//     //         // L_inv_b_view.set(i, val / array[row_diag]);
//     //     }

//     //     // Backward substitution
//     //     //
//     //
//     https://algowiki-project.org/en/Backward_substitution#:~:text=Backward%20substitution%20is%20a%20procedure,is%20a%20lower%20triangular%20matrix.
//     //     for (SizeType i = n - 1; i >= 0; --i) {
//     //         const SizeType row_end = ia[i + 1];
//     //         const SizeType row_diag = idx.idx[i];

//     //         // auto lb_i = lb_view.get(i);
//     //         auto ub_i = ub_view.get(i);

//     //         Scalar val = L_inv_b_view.get(i);

//     //         for (SizeType k = row_diag + 1; k < row_end; ++k) {
//     //             auto j = ja[k];
//     //             val -= array[k] * x_view.get(j);
//     //         }

//     //         // x_view.set(i, std::max(lb_i, std::min(ub_i, val / array[row_diag])));

//     //         x_view.set(i, std::min(ub_i, val / array[row_diag]));
//     //     }
//     // }

//     template <class Block>
//     void block_decompose_aux(const CRSMatrix &mat,
//                              const int block_size,
//                              std::vector<Block> &diag,
//                              CRSMatrix &out,
//                              bool sort_columns = true) {
//         CRSMatrix l_mat;
//         local_block_view(mat, l_mat);

//         out.copy(l_mat);

//         PetscSeqAIJRaw m_raw(out.raw_type());

//         SizeType n = m_raw.n;
//         const SizeType *ia = m_raw.ia;
//         const SizeType *ja = m_raw.ja;
//         Scalar *array = m_raw.array;

//         auto n_blocks = n / block_size;

//         std::vector<SizeType> row_ptr(n_blocks + 1, 0);
//         std::vector<SizeType> block_pattern(n_blocks);

//         for (SizeType block_i = 0; block_i < n_blocks; block_i++) {
//             const SizeType offset_i = block_i * block_size;

//             for (SizeType i = offset_i; i < offset_i + block_size; ++i) {
//                 const SizeType row_begin = ia[i];
//                 const SizeType row_end = ia[i + 1];

//                 for (SizeType k = row_begin; k < row_end; ++k) {
//                     SizeType j = ja[k];
//                     SizeType block_j = j / block_size;

//                     if (block_pattern[block_j] == 0) {
//                         ++row_ptr[block_i + 1];
//                         block_pattern[block_j] = 1.0;
//                     }
//                 }
//             }

//             for (SizeType i = offset_i; i < offset_i + block_size; ++i) {
//                 const SizeType row_begin = ia[i];
//                 const SizeType row_end = ia[i + 1];

//                 for (SizeType k = row_begin; k < row_end; ++k) {
//                     SizeType j = ja[k];
//                     SizeType block_j = j / block_size;
//                     block_pattern[block_j] = 0;
//                 }
//             }
//         }

//         for (SizeType i = 0; i < n_blocks; ++i) {
//             row_ptr[i + 1] += row_ptr[i];
//         }

//         std::vector<SizeType> col_idx(row_ptr[n_blocks], 0);
//         std::vector<SizeType> count(n_blocks, 0);

//         for (SizeType block_i = 0; block_i < n_blocks; block_i++) {
//             const SizeType offset_i = block_i * block_size;

//             for (SizeType i = offset_i; i < offset_i + block_size; ++i) {
//                 const SizeType row_begin = ia[i];
//                 const SizeType row_end = ia[i + 1];

//                 for (SizeType k = row_begin; k < row_end; ++k) {
//                     SizeType j = ja[k];
//                     SizeType block_j = j / block_size;

//                     if (block_pattern[block_j] == 0) {
//                         block_pattern[block_j] = 1.0;
//                         SizeType idx_offset = row_ptr[block_i] + count[block_i];
//                         col_idx[idx_offset] = block_j;
//                         ++count[block_i];
//                     }
//                 }
//             }

//             for (SizeType i = offset_i; i < offset_i + block_size; ++i) {
//                 const SizeType row_begin = ia[i];
//                 const SizeType row_end = ia[i + 1];

//                 for (SizeType k = row_begin; k < row_end; ++k) {
//                     SizeType j = ja[k];
//                     SizeType block_j = j / block_size;
//                     block_pattern[block_j] = 0;
//                 }
//             }

//             if (sort_columns) {
//                 std::sort(&col_idx[row_ptr[block_i]], &col_idx[row_ptr[block_i + 1] - 1]);
//             }
//         }

//         const SizeType block_size_2 = block_size * block_size;
//         std::vector<Scalar> values(row_ptr[n_blocks] * block_size_2, 0.0);

//         for (SizeType block_i = 0; block_i < n_blocks; block_i++) {
//             const SizeType offset_i = block_i * block_size;
//             const SizeType block_row_begin = row_ptr[block_i];
//             const SizeType block_row_end = row_ptr[block_i + 1];

//             for (SizeType i = offset_i; i < offset_i + block_size; ++i) {
//                 const SizeType row_begin = ia[i];
//                 const SizeType row_end = ia[i + 1];
//                 const SizeType sub_i = i - offset_i;

//                 SizeType current_block = block_row_begin;

//                 for (SizeType k = row_begin; k < row_end; ++k) {
//                     const SizeType j = ja[k];
//                     const Scalar value = array[k];

//                     const SizeType block_j = j / block_size;
//                     const SizeType sub_j = j - block_j * block_size;

//                     SizeType k_block = current_block;
//                     for (; k_block < block_row_end; ++k_block) {
//                         if (col_idx[k_block] == block_j) {
//                             current_block = k_block;
//                             break;
//                         }
//                     }

//                     assert(col_idx[k_block] == block_j);

//                     values[k_block * block_size_2 + (sub_i * block_size) + sub_j] = value;
//                 }
//             }
//         }

//         // for (SizeType block_i = 0; block_i < n_blocks; block_i++) {
//         //     const SizeType row_begin = row_ptr[block_i];
//         //     const SizeType row_end = row_ptr[block_i + 1];

//         //     for (SizeType k = row_begin; k < row_end; ++k) {
//         //         std::cout << block_i << " " << col_idx[k] << "\n";

//         //         for (SizeType sub_i = 0; sub_i < block_size; ++sub_i) {
//         //             for (SizeType sub_j = 0; sub_j < block_size; ++sub_j) {
//         //                 std::cout << values[k * block_size_2 + (sub_i * block_size) + sub_j] << " ";
//         //             }

//         //             std::cout << "\n";
//         //         }

//         //         std::cout << "\n";
//         //     }
//         // }
//     }

//     void ILUDecompose<CRSMatrix, HOMEMADE>::block_decompose(const CRSMatrix &mat, CRSMatrix &out, const bool
//     modified) {
//         std::vector<StaticMatrix<Scalar, 2, 2>> diag;
//         block_decompose_aux(mat, 2, diag, out);
//     }

// }  // namespace utopia
