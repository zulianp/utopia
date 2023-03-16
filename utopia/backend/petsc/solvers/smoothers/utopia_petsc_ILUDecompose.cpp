#include "utopia_petsc_ILUDecompose.hpp"
#include "utopia_petsc.hpp"

#include <vector>
#include "utopia_Views.hpp"

#include "utopia_CRSToBlockCRS.hpp"
#include "utopia_ILUDecompose.hpp"

// #include "utopia_petsc_BlockCrsView.hpp"
#include "utopia_petsc_CrsView.hpp"

namespace utopia {

    bool ILUDecompose<PetscMatrix, PETSC>::decompose(const PetscMatrix &mat, PetscMatrix &out, const bool modified) {
        using ScalarView = utopia::ArrayView<PetscScalar>;
        using IndexView = utopia::ArrayView<const PetscInt>;

        PetscMatrix l_mat;
        local_block_view(mat, l_mat);

        // perform copy
        out.copy(l_mat);

        PetscCrsView crs_view(out.raw_type());
        PetscInt n = crs_view.rows();

        auto row_ptr = crs_view.row_ptr();
        auto colidx = crs_view.colidx();
        auto values = crs_view.values();

        CRSMatrix<ScalarView, IndexView> crs(row_ptr, colidx, values, n);

        assert(crs.is_valid());

#ifndef NDEBUG
        const PetscScalar orginal_min = (mat.comm().size() == 1) ? mat.min() : l_mat.min();
        const PetscScalar orginal_max = (mat.comm().size() == 1) ? mat.max() : l_mat.max();

        const PetscScalar crs_min = crs.min();
        const PetscScalar crs_max = crs.max();

        assert(orginal_min == crs_min);
        assert(orginal_max == crs_max);
#endif

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

        PetscCrsView crs_view(ilu.raw_type());
        PetscInt n = crs_view.rows();

        auto row_ptr = crs_view.row_ptr();
        auto colidx = crs_view.colidx();
        auto values = crs_view.values();

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

        PetscCrsView crs_view(l_mat.raw_type());
        PetscInt n = crs_view.rows();

        auto row_ptr = crs_view.row_ptr();
        auto colidx = crs_view.colidx();
        auto values = crs_view.values();

        CRSMatrix<ScalarView, IndexView, 1> crs(row_ptr, colidx, values, n);
        out.row_ptr().clear();
        out.colidx().clear();
        out.values().clear();
        convert(crs, out);

        assert(out.is_valid());

#ifndef NDEBUG
        const PetscScalar orginal_min = (in.comm().size() == 1) ? in.min() : l_mat.min();
        const PetscScalar orginal_max = (in.comm().size() == 1) ? in.max() : l_mat.max();

        const PetscScalar out_min = out.min();
        const PetscScalar out_max = out.max();

        assert(orginal_min == out_min);
        assert(orginal_max == out_max);
#endif
    }

    template <int BlockSize>
    void crs_block_matrix_split_diag(const PetscMatrix &in,
                                     CRSMatrix<std::vector<PetscScalar>, std::vector<PetscInt>, BlockSize> &out,
                                     std::vector<PetscReal> &diag) {
        using ScalarView = utopia::ArrayView<PetscScalar>;
        using IndexView = utopia::ArrayView<const PetscInt>;

        PetscMatrix l_mat;
        local_block_view(in, l_mat);

        PetscCrsView crs_view(l_mat.raw_type());
        PetscInt n = crs_view.rows();
        // PetscInt nnz = crs_view.nnz();

        auto row_ptr = crs_view.row_ptr();
        auto colidx = crs_view.colidx();
        auto values = crs_view.values();

        CRSMatrix<ScalarView, IndexView, 1> crs(row_ptr, colidx, values, n);
        out.row_ptr().clear();
        out.colidx().clear();
        out.values().clear();
        convert_split_diag(crs, out, diag);
    }

    template <int BlockSize>
    void crs_block_matrix_update(const PetscMatrix &in,
                                 CRSMatrix<std::vector<PetscScalar>, std::vector<PetscInt>, BlockSize> &out,
                                 std::vector<PetscReal> &diag) {
        using ScalarView = utopia::ArrayView<PetscScalar>;
        using IndexView = utopia::ArrayView<const PetscInt>;

        PetscMatrix l_mat;
        local_block_view(in, l_mat);

        PetscCrsView crs_view(l_mat.raw_type());
        PetscInt n = crs_view.rows();
        // PetscInt nnz = crs_view.nnz();

        auto row_ptr = crs_view.row_ptr();
        auto colidx = crs_view.colidx();
        auto values = crs_view.values();

        CRSMatrix<ScalarView, IndexView, 1> crs(row_ptr, colidx, values, n);
        convert_split_diag_update(crs, out, diag);
    }
    template void crs_block_matrix<1>(const PetscMatrix &,
                                      CRSMatrix<std::vector<PetscScalar>, std::vector<PetscInt>, 1> &);
    template void crs_block_matrix<2>(const PetscMatrix &,
                                      CRSMatrix<std::vector<PetscScalar>, std::vector<PetscInt>, 2> &);
    template void crs_block_matrix<3>(const PetscMatrix &,
                                      CRSMatrix<std::vector<PetscScalar>, std::vector<PetscInt>, 3> &);
    template void crs_block_matrix<4>(const PetscMatrix &,
                                      CRSMatrix<std::vector<PetscScalar>, std::vector<PetscInt>, 4> &);

    template void crs_block_matrix_split_diag<2>(const PetscMatrix &,
                                                 CRSMatrix<std::vector<PetscScalar>, std::vector<PetscInt>, 2> &,
                                                 std::vector<PetscReal> &);

    template void crs_block_matrix_split_diag<3>(const PetscMatrix &,
                                                 CRSMatrix<std::vector<PetscScalar>, std::vector<PetscInt>, 3> &,
                                                 std::vector<PetscReal> &);

    template void crs_block_matrix_split_diag<4>(const PetscMatrix &,
                                                 CRSMatrix<std::vector<PetscScalar>, std::vector<PetscInt>, 4> &,
                                                 std::vector<PetscReal> &);

    template void crs_block_matrix_update<2>(const PetscMatrix &,
                                             CRSMatrix<std::vector<PetscScalar>, std::vector<PetscInt>, 2> &,
                                             std::vector<PetscReal> &);

    template void crs_block_matrix_update<3>(const PetscMatrix &,
                                             CRSMatrix<std::vector<PetscScalar>, std::vector<PetscInt>, 3> &,
                                             std::vector<PetscReal> &);

    template void crs_block_matrix_update<4>(const PetscMatrix &,
                                             CRSMatrix<std::vector<PetscScalar>, std::vector<PetscInt>, 4> &,
                                             std::vector<PetscReal> &);

    bool ILUDecompose<PetscMatrix, PETSC>::update(const PetscMatrix &mat) { return decompose(mat, ilu_, modified_); }

    void ILUDecompose<PetscMatrix, PETSC>::apply(const PetscVector &b, PetscVector &x) { apply(ilu_, b, x); }

    void ILUDecompose<PetscMatrix, PETSC>::read(Input &in) { in.get("modified", modified_); }

    template <int BlockSize>
    bool BlockILUAlgorithm<PetscMatrix, BlockSize>::update(const PetscMatrix &mat) {
        UTOPIA_TRACE_REGION_BEGIN("BlockILUAlgorithm::update(...)");
        PetscMatrix local_A;
        local_block_view(mat, local_A);
        crs_block_matrix(local_A, ilu_);

        if (print_matrices_) {
            disp(ilu_);
        }

        L_inv_b_.zeros(row_layout(mat));
        bool ok = ilu_decompose(ilu_);

        if (print_matrices_) {
            disp(ilu_);
        }

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
    void BlockILUAlgorithm<PetscMatrix, BlockSize>::read(Input &in) {
        in.get("print_matrices", print_matrices_);
    }

    template class BlockILUAlgorithm<PetscMatrix, 1>;
    template class BlockILUAlgorithm<PetscMatrix, 2>;
    template class BlockILUAlgorithm<PetscMatrix, 3>;
    template class BlockILUAlgorithm<PetscMatrix, 4>;

}  // namespace utopia
