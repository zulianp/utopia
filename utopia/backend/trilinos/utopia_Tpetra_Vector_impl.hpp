#ifndef UTOPIA_TPETRA_VECTOR_IMPL_HPP
#define UTOPIA_TPETRA_VECTOR_IMPL_HPP

#include "utopia_Tpetra_Vector.hpp"

#include <Trilinos_version.h>

#if (TRILINOS_MAJOR_MINOR_VERSION >= 130100 && UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE)
#include <Tpetra_Access.hpp>
#endif

namespace utopia {

    template <class F>
    void TpetraVector::transform_values(F f) {
#if (TRILINOS_MAJOR_MINOR_VERSION >= 130100 && UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE)
        auto k_res = raw_type()->template getLocalView<ExecutionSpace>(Tpetra::Access::ReadWrite);
#else
        auto k_res = raw_type()->template getLocalView<ExecutionSpace>();

#endif

        assert(k_res.extent(0) > 0);

        Kokkos::parallel_for(
            k_res.extent(0), KOKKOS_LAMBDA(const int i) { k_res(i, 0) = f(k_res(i, 0)); });

        Kokkos::fence();
    }

}  // namespace utopia

#endif  // UTOPIA_TPETRA_VECTOR_IMPL_HPP
