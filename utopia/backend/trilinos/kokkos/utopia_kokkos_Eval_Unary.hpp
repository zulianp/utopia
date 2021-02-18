#ifndef UTOPIA_KOKKOS_EVAL_UNARY_HPP
#define UTOPIA_KOKKOS_EVAL_UNARY_HPP

#include <Kokkos_Core.hpp>
#include <cassert>
#include "utopia_kokkos_Operations.hpp"

namespace utopia {

    template <class Vector, class Op>
    class KokkosEvalUnary {
    public:
        inline static void eval(const Op op, Vector &vec) {
            using ExecutionSpaceT = typename Vector::ExecutionSpace;
            using Scalar = typename Vector::Scalar;

            auto k_res = vec.implementation().template getLocalView<ExecutionSpaceT>();

            assert(k_res.extent(0) > 0);

            KokkosOp<Scalar, Op> k_op(op);

            Kokkos::parallel_for(
                k_res.extent(0), KOKKOS_LAMBDA(const int i) { k_res(i, 0) = k_op.apply(k_res(i, 0)); });
        }
    };

}  // namespace utopia

#endif  // UTOPIA_KOKKOS_EVAL_UNARY_HPP
