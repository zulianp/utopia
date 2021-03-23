#include "utopia_moonolith_ConvertTensor.hpp"

#include "moonolith_sparse_matrix.hpp"

namespace utopia {

    void ConvertTensor<::moonolith::SparseMatrix<UScalar>, USparseMatrix, 2>::apply(
        const ::moonolith::SparseMatrix<UScalar> &in,
        USparseMatrix &out) {
        USizeType n_local_rows = in.local_rows();
        USizeType n_local_cols = in.local_cols();

        USizeType rows = in.rows();
        USizeType cols = in.cols();

        auto &&comm = out.comm();

        USizeType begin = 0;
        comm.exscan_sum(&n_local_rows, &begin, 1);

        UIndexArray d_nnz(n_local_rows, 0), o_nnz(n_local_rows, 0);

        assert(comm.rank() != 0 || begin == 0);

        Range r(begin, begin + n_local_rows);

        for (auto it = in.iter(); it; ++it) {
            const UScalar row = it.row() - begin;
            if (r.inside(it.col())) {
                ++d_nnz[row];
            } else {
                ++o_nnz[row];
            }
        }

        auto lo = layout(comm, n_local_rows, n_local_cols, rows, cols);
        out.sparse(lo, d_nnz, o_nnz);

        {
            Write<USparseMatrix> write(out);
            for (auto it = in.iter(); it; ++it) {
                out.set(it.row(), it.col(), *it);
            }
        }
    }

    void ConvertTensor<::moonolith::SparseMatrix<UScalar>, UVector, 1>::apply(
        const ::moonolith::SparseMatrix<UScalar> &in,
        UVector &out) {
        auto n_local_rows = in.local_rows();
        auto n_local_cols = in.local_cols();

        assert(n_local_cols == 1);

        out.zeros(layout(in.comm().get(), n_local_rows, Traits<UVector>::determine()));

        {
            Write<UVector> write(out);
            for (auto it = in.iter(); it; ++it) {
                assert(it.col() == 0);
                out.set(it.row(), *it);
            }
        }
    }

}  // namespace utopia
