#include "utopia_petsc_ILUDecompose.hpp"
#include "utopia_petsc.hpp"

#include <vector>
#include "utopia_Views.hpp"

#include "utopia_CRSToBlockCRS.hpp"
#include "utopia_ILUDecompose.hpp"

namespace utopia {

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

    bool ILUDecompose<PetscMatrix, PETSC>::decompose(const PetscMatrix &mat, PetscMatrix &out, const bool modified) {
        using ScalarView = utopia::ArrayView<PetscScalar>;
        using IndexView = utopia::ArrayView<const PetscInt>;

        PetscMatrix l_mat;
        local_block_view(mat, l_mat);

        // perform copy
        out.copy(l_mat);

        PetscSeqAIJRaw m_raw(out.raw_type());
        PetscInt n = m_raw.n;
        PetscInt nnz = m_raw.ia[n];

        const PetscInt *ia = m_raw.ia;
        const PetscInt *ja = m_raw.ja;
        PetscScalar *array = m_raw.array;

        ScalarView values(array, nnz);
        IndexView row_ptr(ia, n + 1);
        IndexView colidx(ja, nnz);

        CRSMatrix<ScalarView, IndexView> crs(row_ptr, colidx, values, n);
        return ilu_decompose(crs, modified);
    }

    void ILUDecompose<PetscMatrix, PETSC>::apply(const PetscMatrix &ilu, const PetscVector &b, PetscVector &x) {
        using ScalarView = utopia::ArrayView<PetscScalar>;
        using IndexView = utopia::ArrayView<const PetscInt>;

        PetscVector L_inv_b(layout(b), 0.0);
        // x.set(0.0);

        auto b_view = const_local_view_device(b);
        auto L_inv_b_view = local_view_device(L_inv_b);
        auto x_view = local_view_device(x);

        PetscSeqAIJRaw m_raw(ilu.raw_type());

        PetscInt n = m_raw.n;
        PetscInt nnz = m_raw.ia[n];

        const PetscInt *ia = m_raw.ia;
        const PetscInt *ja = m_raw.ja;
        PetscScalar *array = m_raw.array;

        ScalarView values(array, nnz);
        IndexView row_ptr(ia, n + 1);
        IndexView colidx(ja, nnz);

        CRSMatrix<ScalarView, IndexView> crs(row_ptr, colidx, values, n);

        auto b_array = b_view.array();
        auto L_inv_b_array = L_inv_b_view.array();
        auto x_array = x_view.array();

        ilu_apply(crs, b_array, L_inv_b_array, x_array);
    }

    template <int BlockSize>
    void crs_block_matrix(const PetscMatrix &in,
                          CRSMatrix<std::vector<PetscScalar>, std::vector<PetscInt>, BlockSize> &out) {
        using ScalarView = utopia::ArrayView<PetscScalar>;
        using IndexView = utopia::ArrayView<const PetscInt>;

        PetscMatrix l_mat;
        local_block_view(in, l_mat);

        // perform copy
        // out.copy(l_mat);

        PetscSeqAIJRaw m_raw(l_mat.raw_type());
        PetscInt n = m_raw.n;
        PetscInt nnz = m_raw.ia[n];

        const PetscInt *ia = m_raw.ia;
        const PetscInt *ja = m_raw.ja;
        PetscScalar *array = m_raw.array;

        ScalarView values(array, nnz);
        IndexView row_ptr(ia, n + 1);
        IndexView colidx(ja, nnz);

        CRSMatrix<ScalarView, IndexView, 1> crs(row_ptr, colidx, values, n);
        out.row_ptr().clear();
        out.colidx().clear();
        out.values().clear();
        convert(crs, out);
    }

    template void crs_block_matrix<2>(const PetscMatrix &,
                                      CRSMatrix<std::vector<PetscScalar>, std::vector<PetscInt>, 2> &);
    template void crs_block_matrix<3>(const PetscMatrix &,
                                      CRSMatrix<std::vector<PetscScalar>, std::vector<PetscInt>, 3> &);
    template void crs_block_matrix<4>(const PetscMatrix &,
                                      CRSMatrix<std::vector<PetscScalar>, std::vector<PetscInt>, 4> &);

    bool ILUDecompose<PetscMatrix, PETSC>::update(const PetscMatrix &mat) { return decompose(mat, ilu_, modified_); }

    void ILUDecompose<PetscMatrix, PETSC>::apply(const PetscVector &b, PetscVector &x) { apply(ilu_, b, x); }

    void ILUDecompose<PetscMatrix, PETSC>::read(Input &in) { in.get("modified", modified_); }

    template <int BlockSize>
    bool BlockILUAlgorithm<PetscMatrix, BlockSize>::update(const PetscMatrix &mat) {
        UTOPIA_TRACE_REGION_BEGIN("BlockILUAlgorithm::update(...)");
        PetscMatrix local_A;
        local_block_view(mat, local_A);
        crs_block_matrix(local_A, ilu_);

        L_inv_b_.zeros(row_layout(mat));
        bool ok = ilu_decompose(ilu_);

        UTOPIA_TRACE_REGION_END("BlockILUAlgorithm::update(...)");
        return ok;
    }

    template <int BlockSize>
    void BlockILUAlgorithm<PetscMatrix, BlockSize>::apply(const PetscVector &b, PetscVector &x) {
        UTOPIA_TRACE_REGION_BEGIN("BlockILUAlgorithm::apply(...)");
        auto b_view = const_local_view_device(b);
        auto L_inv_b_view = local_view_device(L_inv_b_);
        auto x_view = local_view_device(x);

        auto b_array = b_view.array();
        auto x_array = x_view.array();
        auto L_inv_b_array = L_inv_b_view.array();

        ilu_apply(ilu_, b_array, L_inv_b_array, x_array);

        UTOPIA_TRACE_REGION_END("BlockILUAlgorithm::apply(...)");
    }

    template <int BlockSize>
    void BlockILUAlgorithm<PetscMatrix, BlockSize>::read(Input &) {}

    template class BlockILUAlgorithm<PetscMatrix, 2>;
    template class BlockILUAlgorithm<PetscMatrix, 3>;
    template class BlockILUAlgorithm<PetscMatrix, 4>;

}  // namespace utopia
