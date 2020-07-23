#include "utopia_Base.hpp"

#include "utopia_kokkos_Eval_Reduce.hpp"
#include "utopia_kokkos_Operations.hpp"
#include "utopia_trilinos_MaxRowNNZ.hpp"

#include <Kokkos_View.hpp>

namespace utopia {

    Traits<TpetraMatrix>::SizeType MaxRowNNZ<TpetraMatrix, TRILINOS>::apply(const TpetraMatrix &in) {
        using ExecutionSpaceT = typename TpetraVector::ExecutionSpace;
        using LocalMatrix = TpetraMatrix::CrsMatrixType::local_matrix_type;
        using ViewType = Kokkos::View<LocalSizeType *[1], ExecutionSpaceT>;

        auto r = in.row_range();
        if (r.empty()) {
            return 0;
        }

        auto impl = in.raw_type();
        auto local_mat = impl->getLocalMatrix();

        auto n = local_mat.numRows();

        auto row_map = impl->getRowMap()->getLocalMap();
        auto col_map = impl->getColMap()->getLocalMap();

        ViewType nnz("nnz", r.extent());

        Kokkos::parallel_for(
            "MaxRowNNZ::apply", n, UTOPIA_LAMBDA(const int &i) { nnz(i, 0) = local_mat.row(i).length; });

        Scalar ret = 0;
        KokkosOp<Scalar, Max> kop;
        OpFunctor<ViewType, KokkosOp<Scalar, Max>, Scalar> functor{kop, nnz, ret};
        Kokkos::parallel_reduce(nnz.extent(0), functor, ret);
        return ret;
    }
}  // namespace utopia