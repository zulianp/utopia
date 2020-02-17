#ifndef UTOPIA_KOKKOS_EVAL_BINARY_HPP
#define UTOPIA_KOKKOS_EVAL_BINARY_HPP

#include "utopia_kokkos_Operations.hpp"
#include <Kokkos_Core.hpp>
#include <iostream>
#include <cassert>

namespace utopia {
    template<class Vector, class Op>
    class KokkosEvalBinary {
    public:
        using Scalar = typename Traits<Vector>::Scalar;

        inline static void eval(const Vector &lhs, const Op op, const Vector &rhs, Vector &result)
        {
            using ExecutionSpaceT = typename Vector::ExecutionSpace;

            assert(!lhs.empty());
            assert(!rhs.empty());
            assert(rhs.size() == lhs.size());
            assert(rhs.local_size() == lhs.local_size());

            if(result.empty() || result.size() != rhs.size())
            {
                result.init(rhs.implementation().getMap());
            }

            auto k_lhs = lhs.implementation().template getLocalView<ExecutionSpaceT> ();
            auto k_rhs = rhs.implementation().template getLocalView<ExecutionSpaceT> ();
            auto k_res = result.implementation().template getLocalView<ExecutionSpaceT> ();

            KokkosOp<Scalar, Op> k_op;
            Kokkos::parallel_for(k_lhs.extent(0), KOKKOS_LAMBDA (const int i) {
                k_res(i, 0) = k_op.apply(k_lhs(i, 0), k_rhs(i, 0));
            });
        }

        inline static void eval(const Vector &lhs, const Op op, const Scalar &rhs, Vector &result)
        {
            using ExecutionSpaceT = typename Vector::ExecutionSpace;

            assert(!lhs.empty());

            if(result.empty())
            {
                result.init(lhs.implementation().getMap());
            }

            auto k_lhs = lhs.implementation().template getLocalView<ExecutionSpaceT> ();
            auto k_res = result.implementation().template getLocalView<ExecutionSpaceT> ();

            KokkosOp<Scalar, Op> k_op;
            Kokkos::parallel_for(k_lhs.extent(0), KOKKOS_LAMBDA (const int i) {
                k_res(i, 0) = k_op.apply(k_lhs(i, 0), rhs);
            });
        }
    };
}

#endif

