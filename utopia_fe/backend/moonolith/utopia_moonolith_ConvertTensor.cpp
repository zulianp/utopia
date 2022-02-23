#include "utopia_moonolith_ConvertTensor.hpp"

#include "moonolith_sparse_matrix.hpp"

namespace utopia {

    void ConvertTensor<::moonolith::SparseMatrix<UScalar>, USparseMatrix, 2>::apply(
        const ::moonolith::SparseMatrix<UScalar> &in,
        USparseMatrix &out) {
        UTOPIA_TRACE_REGION_BEGIN("ConvertTensor::apply | allocate crs");

        USizeType n_local_rows = in.local_rows();
        USizeType n_local_cols = in.local_cols();

        USizeType rows = in.rows();
        USizeType cols = in.cols();

        assert(n_local_rows >= 0);
        assert(n_local_cols >= 0);
        assert(rows > 0);
        assert(cols > 0);

        auto &&comm = out.comm();

        USizeType local_offsets[2] = {n_local_rows, n_local_cols};
        USizeType begins[2] = {0, 0};
        comm.exscan_sum(local_offsets, begins, 2);

        UIndexArray d_nnz(n_local_rows, 0), o_nnz(n_local_rows, 0);

        assert(comm.rank() != 0 || begins[0] == 0);
        assert(comm.rank() != 0 || begins[1] == 0);

        Range r(begins[0], begins[0] + n_local_rows);
        Range cr(begins[1], begins[1] + n_local_cols);

        for (auto it = in.iter(); it; ++it) {
            const UScalar row = it.row() - begins[0];
            assert(r.inside(it.row()));

            if (cr.inside(it.col())) {
                ++d_nnz[row];
            } else {
                ++o_nnz[row];
            }
        }

        auto lo = layout(comm, n_local_rows, n_local_cols, rows, cols);
        out.sparse(lo, d_nnz, o_nnz);

        UTOPIA_TRACE_REGION_END("ConvertTensor::apply | allocate crs");

        UTOPIA_TRACE_REGION_BEGIN("ConvertTensor::apply | fill crs");

        {
            Write<USparseMatrix> write(out);
            for (auto it = in.iter(); it; ++it) {
                out.set(it.row(), it.col(), *it);
            }
        }

        UTOPIA_TRACE_REGION_END("ConvertTensor::apply | fill crs");
    }

    void ConvertTensor<::moonolith::SparseMatrix<UScalar>, UVector, 1>::apply(
        const ::moonolith::SparseMatrix<UScalar> &in,
        UVector &out) {
        auto n_local_rows = in.local_rows();

        assert(in.local_cols() == 1);
        assert(n_local_rows >= 0);
        assert(in.rows() > 0);

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
