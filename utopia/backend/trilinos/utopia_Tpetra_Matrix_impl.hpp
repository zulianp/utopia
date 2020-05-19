#ifndef UTOPIA_TPETRA_MATRIX_IMPL_HPP
#define UTOPIA_TPETRA_MATRIX_IMPL_HPP

#include "utopia_Tpetra_Matrix.hpp"
#include "utopia_kokkos_Eval_Reduce.hpp"

namespace utopia {

    template <class Op>
    TpetraMatrix::Scalar TpetraMatrix::parallel_reduce_values(Op op, const Scalar &initial_value) const {
        using LocalMatrix = typename CrsMatrixType::local_matrix_type;
        using Data = typename LocalMatrix::values_type;
        const LocalMatrix &local_mat = raw_type()->getLocalMatrix();
        const Data &data = local_mat.values;

        Scalar ret = initial_value;
        KokkosOp<Scalar, Op> kop;
        OpFunctor<Data, KokkosOp<Scalar, Op>, Scalar> functor{kop, data, initial_value};
        Kokkos::parallel_reduce(data.extent(0), functor, ret);

        ret = comm().reduce(op, ret);
        return ret;
    }
}  // namespace utopia

#endif  // UTOPIA_TPETRA_MATRIX_IMPL_HPP
