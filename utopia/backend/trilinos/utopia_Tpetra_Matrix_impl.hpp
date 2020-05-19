#ifndef UTOPIA_TPETRA_MATRIX_IMPL_HPP
#define UTOPIA_TPETRA_MATRIX_IMPL_HPP

#include "utopia_Tpetra_Matrix.hpp"
#include "utopia_kokkos_Eval_Reduce.hpp"

namespace utopia {

    template <class Op>
    TpetraMatrix::Scalar TpetraMatrix::local_parallel_reduce_values(Op op, const Scalar &initial_value) const {
        using LocalMatrix = typename CrsMatrixType::local_matrix_type;
        using Data = typename LocalMatrix::values_type;
        const LocalMatrix &local_mat = raw_type()->getLocalMatrix();
        const Data &data = local_mat.values;

        Scalar ret = initial_value;
        KokkosOp<Scalar, Op> kop;
        OpFunctor<Data, KokkosOp<Scalar, Op>, Scalar> functor{kop, data, initial_value};
        Kokkos::parallel_reduce(data.extent(0), functor, ret);
        return ret;
    }

    template <class Op, class MPIOp>
    TpetraMatrix::Scalar TpetraMatrix::parallel_reduce_values(Op op, MPIOp mpi_op, const Scalar &initial_value) const {
        auto ret = local_parallel_reduce_values(op, initial_value);
        ret = comm().reduce(mpi_op, ret);
        return ret;
    }

    template <class Op>
    void TpetraMatrix::aux_transform(const Op &op) {
        KokkosOp<Scalar, Op> k_op(op);
        transform_values(KOKKOS_LAMBDA(const Scalar value)->Scalar { return k_op.apply(value); });
    }

    template <class F>
    void TpetraMatrix::transform_values(F op) {
        using LocalMatrix = typename CrsMatrixType::local_matrix_type;
        using Data = typename LocalMatrix::values_type;
        const LocalMatrix &local_mat = raw_type()->getLocalMatrix();
        const Data &data = local_mat.values;

        Kokkos::parallel_for(data.extent(0), KOKKOS_LAMBDA(const SizeType i) { data(i) = op(data(i)); });
    }

    template <class Op>
    void TpetraMatrix::transform_ijv(Op op) {
        typedef Kokkos::TeamPolicy<> team_policy;
        typedef Kokkos::TeamPolicy<>::member_type member_type;

        auto r = this->row_range();

        if (r.empty()) {
            return;
        }

        auto impl = this->raw_type();
        auto local_mat = impl->getLocalMatrix();

        auto n = local_mat.numRows();

        auto row_map = impl->getRowMap()->getLocalMap();
        auto col_map = impl->getColMap()->getLocalMap();

        Kokkos::parallel_for(
            "TpetraMatrix::transform_ijv", team_policy(n, Kokkos::AUTO), KOKKOS_LAMBDA(const member_type &team_member) {
                const int row_ind = team_member.league_rank();
                auto row = local_mat.row(row_ind);
                auto n_values = row.length;

                Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, n_values), [&](const int k) {
                    auto &val = row.value(k);
                    auto col_ind = row.colidx(k);
                    val = op(row_map.getGlobalElement(row_ind), col_map.getGlobalElement(col_ind), val);
                });
            });
    }

}  // namespace utopia

#endif  // UTOPIA_TPETRA_MATRIX_IMPL_HPP
