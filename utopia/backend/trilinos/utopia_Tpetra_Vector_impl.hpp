#ifndef UTOPIA_TPETRA_VECTOR_IMPL_HPP
#define UTOPIA_TPETRA_VECTOR_IMPL_HPP

#include "utopia_Tpetra_Vector.hpp"

namespace utopia {

    template <class F>
    void TpetraVector::transform_values(F f) {
        auto k_res = raw_type()->template getLocalView<ExecutionSpace>();

        assert(k_res.extent(0) > 0);

        Kokkos::parallel_for(
            k_res.extent(0), KOKKOS_LAMBDA(const int i) { k_res(i, 0) = f(k_res(i, 0)); });
    }

}  // namespace utopia

#endif  // UTOPIA_TPETRA_VECTOR_IMPL_HPP
