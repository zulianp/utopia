#include "utopia_moonolith_ConvertTensor.hpp"

#include "moonolith_sparse_matrix.hpp"

namespace utopia {

    void ConvertTensor<::moonolith::SparseMatrix<UScalar>, USparseMatrix, 2>::apply(
        const ::moonolith::SparseMatrix<UScalar> &in,
        USparseMatrix &out) {}

    void ConvertTensor<::moonolith::SparseMatrix<UScalar>, UVector, 1>::apply(
        const ::moonolith::SparseMatrix<UScalar> &in,
        UVector &out) {}

}  // namespace utopia