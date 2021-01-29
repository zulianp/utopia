#include "utopia_petsc_ILUDecompose.hpp"
#include "utopia_petsc.hpp"

#include <vector>

namespace utopia {

    class DiagIdx {
    public:
        void init(PetscInt n, const PetscInt *ia, const PetscInt *ja) {
            idx.resize(n);
            this->ia = ia;
            this->ja = ja;

            for (PetscInt r = 0; r < n - 1; ++r) {
                const PetscInt row_end = ia[r + 1];

                PetscInt k;
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

        auto find_before_diag(const PetscInt i, const PetscInt j) const -> PetscInt {
            PetscInt search_end = idx[i];
            PetscInt row_begin = ia[i];

            // maybe binary search?
            for (PetscInt k = row_begin; k < search_end; ++k) {
                if (ja[k] == j) {
                    return k;
                }
            }

            return -1;
        };

        auto find_after_diag(const PetscInt i, const PetscInt j) const -> PetscInt {
            PetscInt search_begin = idx[i];
            PetscInt row_end = ia[i + 1];

            // maybe binary search?
            for (PetscInt k = search_begin; k < row_end; ++k) {
                if (ja[k] == j) {
                    return k;
                }
            }

            return -1;
        };

        auto find(const PetscInt i, const PetscInt j) const -> PetscInt {
            if (i == j) return idx[i];

            if (j < i) {
                return find_before_diag(i, j);
            } else {
                return find_after_diag(i, j);
            }
        };

        std::vector<PetscInt> idx;
        const PetscInt *ia;
        const PetscInt *ja;
    };

    class PetscSeqAIJRaw {
    public:
        PetscSeqAIJRaw(Mat raw_mat) : raw_mat(raw_mat) {
            err = MatGetRowIJ(raw_mat, 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, &ja, &done);
            assert(err == 0);
            assert(done == PETSC_TRUE);

            if (!done) {
                utopia::err()
                    << "PetscMatrix::read_petsc_seqaij_impl(const Op &op): MatGetRowIJ failed to provide what "
                       "was asked.\n";
                abort();
            }

            MatSeqAIJGetArray(raw_mat, &array);
        }

        ~PetscSeqAIJRaw() {
            MatSeqAIJRestoreArray(raw_mat, &array);
            err = MatRestoreRowIJ(raw_mat, 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, &ja, &done);
            assert(err == 0);
            UTOPIA_UNUSED(err);
        }

        Mat raw_mat;
        PetscInt n = 0;
        const PetscInt *ia{nullptr};
        const PetscInt *ja{nullptr};
        PetscScalar *array{nullptr};
        PetscBool done;
        PetscErrorCode err{0};
    };

    template <typename Apply_IJ>
    void decompose_aux(Apply_IJ apply_ij, const PetscMatrix &mat, PetscMatrix &out) {
        PetscMatrix l_mat;
        local_block_view(mat, l_mat);

        // perform copy
        out.copy(l_mat);

        PetscSeqAIJRaw m_raw(out.raw_type());
        PetscInt n = m_raw.n;
        const PetscInt *ia = m_raw.ia;
        const PetscInt *ja = m_raw.ja;
        PetscScalar *array = m_raw.array;

        DiagIdx idx;
        idx.init(n, ia, ja);

        for (PetscInt r = 0; r < n - 1; ++r) {
            const PetscInt row_end = ia[r + 1];
            const PetscInt k = idx.idx[r];
            const PetscScalar d = 1. / array[k];

            for (PetscInt ki = k + 1; ki < row_end; ++ki) {
                auto i = ja[ki];

                auto ir = idx.find_before_diag(i, r);
                if (ir == -1) continue;

                auto ii = idx.idx[i];

                PetscScalar e = (array[ir] *= d);

                for (PetscInt rj = k + 1; rj < row_end; ++rj) {
                    auto j = ja[rj];

                    auto ij = idx.find(i, j);
                    apply_ij(ii, ij, rj, e, array);
                }
            }
        }
    }

    void ILUDecompose<PetscMatrix, PETSC>::decompose(const PetscMatrix &mat, PetscMatrix &out, const bool modified) {
        if (modified) {
            decompose_aux(
                [](const PetscInt ii, const PetscInt ij, const PetscInt rj, const PetscScalar &e, PetscScalar *a) {
                    if (ij != -1) {
                        a[ij] -= e * a[rj];
                    } else {
                        a[ii] -= e * a[rj];
                    }
                },
                mat,
                out);
        } else {
            decompose_aux(
                [](const PetscInt, const PetscInt ij, const PetscInt rj, const PetscScalar &e, PetscScalar *a) {
                    if (ij != -1) {
                        a[ij] -= e * a[rj];
                    }
                },
                mat,
                out);
        }
    }

    void ILUDecompose<PetscMatrix, PETSC>::apply(const PetscMatrix &ilu, const PetscVector &b, PetscVector &x) {
        PetscVector L_inv_b(layout(b), 0.0);
        x.set(0.0);

        auto b_view = const_local_view_device(b);
        auto L_inv_b_view = local_view_device(L_inv_b);
        auto x_view = local_view_device(x);

        PetscSeqAIJRaw m_raw(ilu.raw_type());

        PetscInt n = m_raw.n;
        const PetscInt *ia = m_raw.ia;
        const PetscInt *ja = m_raw.ja;
        PetscScalar *array = m_raw.array;

        DiagIdx idx;
        idx.init(n, ia, ja);

        // Forward substitution
        // https://algowiki-project.org/en/Forward_substitution#:~:text=Forward%20substitution%20is%20the%20process,math%5DL%5B%2Fmath%5D.
        for (PetscInt i = 0; i < n; ++i) {
            const PetscInt row_begin = ia[i];
            const PetscInt row_diag = idx.idx[i];

            PetscScalar val = b_view.get(i);

            for (PetscInt k = row_begin; k < row_diag; ++k) {
                auto j = ja[k];
                val -= array[k] * L_inv_b_view.get(j);
            }

            L_inv_b_view.set(i, val / array[row_diag]);
        }

        // Backward substitution
        // https://algowiki-project.org/en/Backward_substitution#:~:text=Backward%20substitution%20is%20a%20procedure,is%20a%20lower%20triangular%20matrix.
        for (PetscInt i = n - 1; i >= 0; --i) {
            const PetscInt row_end = ia[i + 1];
            const PetscInt row_diag = idx.idx[i];

            PetscScalar val = L_inv_b_view.get(i);

            for (PetscInt k = row_diag + 1; k < row_end; ++k) {
                auto j = ja[k];
                val -= array[k] * x_view.get(j);
            }

            x_view.set(i, val / array[row_diag]);
        }
    }

    // void ILUDecompose<PetscMatrix, PETSC>::apply_vi(const PetscMatrix &ilu,
    //                                                 const PetscVector &lb,
    //                                                 const PetscVector &ub,
    //                                                 const PetscVector &b,
    //                                                 PetscVector &x) {
    //     PetscVector L_inv_b(layout(b), 0.0), U_ub(layout(b), 0.0);
    //     x.set(0.0);

    //     auto b_view = const_local_view_device(b);
    //     auto lb_view = const_local_view_device(lb);
    //     auto ub_view = const_local_view_device(ub);

    //     auto U_ub_view = local_view_device(U_ub);
    //     auto L_inv_b_view = local_view_device(L_inv_b);
    //     auto x_view = local_view_device(x);

    //     PetscSeqAIJRaw m_raw(ilu.raw_type());

    //     PetscInt n = m_raw.n;
    //     const PetscInt *ia = m_raw.ia;
    //     const PetscInt *ja = m_raw.ja;
    //     PetscScalar *array = m_raw.array;

    //     DiagIdx idx;
    //     idx.init(n, ia, ja);

    //     // Backward substitution
    //     for (PetscInt i = n - 1; i >= 0; --i) {
    //         const PetscInt row_end = ia[i + 1];
    //         const PetscInt row_diag = idx.idx[i];

    //         PetscScalar val = ub_view.get(i);

    //         for (PetscInt k = row_diag + 1; k < row_end; ++k) {
    //             auto j = ja[k];
    //             val -= array[k] * U_ub_view.get(j);
    //         }

    //         U_ub_view.set(i, val / array[row_diag]);
    //     }

    //     // Forward substitution
    //     //
    //     https://algowiki-project.org/en/Forward_substitution#:~:text=Forward%20substitution%20is%20the%20process,math%5DL%5B%2Fmath%5D.
    //     for (PetscInt i = 0; i < n; ++i) {
    //         const PetscInt row_begin = ia[i];
    //         const PetscInt row_diag = idx.idx[i];

    //         // auto lb_i = lb_view.get(i);
    //         auto ub_i = U_ub_view.get(i);

    //         PetscScalar val = b_view.get(i);

    //         for (PetscInt k = row_begin; k < row_diag; ++k) {
    //             auto j = ja[k];
    //             val -= array[k] * L_inv_b_view.get(j);
    //         }

    //         // L_inv_b_view.set(i, std::max(lb_i, std::min(ub_i, val / array[row_diag])));

    //         L_inv_b_view.set(i, std::min(ub_i, val / array[row_diag]));

    //         // L_inv_b_view.set(i, val / array[row_diag]);
    //     }

    //     // Backward substitution
    //     //
    //     https://algowiki-project.org/en/Backward_substitution#:~:text=Backward%20substitution%20is%20a%20procedure,is%20a%20lower%20triangular%20matrix.
    //     for (PetscInt i = n - 1; i >= 0; --i) {
    //         const PetscInt row_end = ia[i + 1];
    //         const PetscInt row_diag = idx.idx[i];

    //         // auto lb_i = lb_view.get(i);
    //         auto ub_i = ub_view.get(i);

    //         PetscScalar val = L_inv_b_view.get(i);

    //         for (PetscInt k = row_diag + 1; k < row_end; ++k) {
    //             auto j = ja[k];
    //             val -= array[k] * x_view.get(j);
    //         }

    //         // x_view.set(i, std::max(lb_i, std::min(ub_i, val / array[row_diag])));

    //         x_view.set(i, std::min(ub_i, val / array[row_diag]));
    //     }
    // }

}  // namespace utopia
