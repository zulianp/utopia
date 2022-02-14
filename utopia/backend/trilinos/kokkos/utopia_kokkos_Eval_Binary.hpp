#ifndef UTOPIA_KOKKOS_EVAL_BINARY_HPP
#define UTOPIA_KOKKOS_EVAL_BINARY_HPP

#include <Kokkos_Core.hpp>
#include <cassert>
#include <iostream>
#include "utopia_kokkos_Operations.hpp"

#include <Trilinos_version.h>

#if TRILINOS_MAJOR_VERSION >= 13
#include <Tpetra_Access.hpp>
#endif

namespace utopia {
    template <class Vector, class Op>
    class KokkosEvalBinary {
    public:
        using Scalar = typename Traits<Vector>::Scalar;

        inline static void eval(const Vector &lhs, const Op &, const Vector &rhs, Vector &result) {
            using ExecutionSpaceT = typename Vector::ExecutionSpace;

            assert(!lhs.empty());
            assert(!rhs.empty());
            assert(rhs.size() == lhs.size());
            assert(rhs.local_size() == lhs.local_size());

            if (result.empty() || result.size() != rhs.size()) {
                result.init(rhs.implementation().getMap());
            }

#if TRILINOS_MAJOR_VERSION >= 13
            auto k_lhs = lhs.implementation().template getLocalView<ExecutionSpaceT>(Tpetra::Access::ReadOnly);
            auto k_rhs = rhs.implementation().template getLocalView<ExecutionSpaceT>(Tpetra::Access::ReadOnly);
            auto k_res = result.implementation().template getLocalView<ExecutionSpaceT>(Tpetra::Access::ReadWrite);
#else
            auto k_lhs = lhs.implementation().template getLocalView<ExecutionSpaceT>();
            auto k_rhs = rhs.implementation().template getLocalView<ExecutionSpaceT>();
            auto k_res = result.implementation().template getLocalView<ExecutionSpaceT>();
#endif

            KokkosOp<Scalar, Op> k_op;
            Kokkos::parallel_for(
                k_lhs.extent(0), KOKKOS_LAMBDA(const int i) { k_res(i, 0) = k_op.apply(k_lhs(i, 0), k_rhs(i, 0)); });

            Kokkos::fence();
        }

        inline static void eval(const Vector &lhs, const Op &, const Scalar &rhs, Vector &result) {
            using ExecutionSpaceT = typename Vector::ExecutionSpace;

            assert(!lhs.empty());

            if (result.empty()) {
                result.init(lhs.implementation().getMap());
            }

#if TRILINOS_MAJOR_VERSION >= 13
            auto k_lhs = lhs.implementation().template getLocalView<ExecutionSpaceT>(Tpetra::Access::ReadOnly);
            auto k_res = result.implementation().template getLocalView<ExecutionSpaceT>(Tpetra::Access::ReadWrite);
#else
            auto k_lhs = lhs.implementation().template getLocalView<ExecutionSpaceT>();
            auto k_res = result.implementation().template getLocalView<ExecutionSpaceT>();
#endif

            KokkosOp<Scalar, Op> k_op;
            Kokkos::parallel_for(
                k_lhs.extent(0), KOKKOS_LAMBDA(const int i) { k_res(i, 0) = k_op.apply(k_lhs(i, 0), rhs); });

            Kokkos::fence();
        }
    };
}  // namespace utopia

#endif
