#ifndef UTOPIA_KOKKOS_EVAL_REDUCE_HPP
#define UTOPIA_KOKKOS_EVAL_REDUCE_HPP

#include "utopia_Traits.hpp"
#include "utopia_kokkos_Operations.hpp"

#include <Kokkos_Core.hpp>

#include <cassert>
#include <iostream>

#include <Trilinos_version.h>
#include <Tpetra_Access.hpp>

namespace utopia {
    template <typename Data, typename KokkosOp, typename Scalar>
    struct OpFunctor {
        KOKKOS_INLINE_FUNCTION void join(volatile Scalar &val, const volatile Scalar &other) const {
            // Kokkos forces us to have the input values being declared volatile. Hence we need to make copies for the
            // reduction operations
            const Scalar tmp1 = val, tmp2 = other;
            val = op_.apply(tmp1, tmp2);
        }

        KOKKOS_INLINE_FUNCTION void operator()(const int &i, Scalar &val) const { val = op_.apply(val, data_(i, 0)); }

        KOKKOS_INLINE_FUNCTION void init(Scalar &val) const { val = initial_value_; }

        const KokkosOp op_;
        const Data data_;
        const Scalar initial_value_;
    };

    template <class Vector, class Op>
    class KokkosEvalReduce {
    public:
        using Scalar = typename Traits<Vector>::Scalar;

        inline static Scalar eval(const Vector &vec, const Op &, const Scalar &initial_value) {
            using ExecutionSpaceT = typename Vector::ExecutionSpace;
            using Scalar = typename Vector::Scalar;

            assert(!vec.empty());

#if TRILINOS_MAJOR_VERSION >= 13
            using Data = decltype(vec.raw_type()->template getLocalView<ExecutionSpaceT>(Tpetra::Access::ReadOnly));
            auto data = vec.raw_type()->template getLocalView<ExecutionSpaceT>(Tpetra::Access::ReadOnly);
#else
            using Data = decltype(vec.raw_type()->template getLocalView<ExecutionSpaceT>());
            auto data = vec.raw_type()->template getLocalView<ExecutionSpaceT>();
#endif

            Scalar ret = initial_value;
            KokkosOp<Scalar, Op> kop;
            OpFunctor<Data, KokkosOp<Scalar, Op>, Scalar> functor{kop, data, initial_value};
            Kokkos::parallel_reduce(data.extent(0), functor, ret);
            Kokkos::fence();

            return ret;
        }
    };

}  // namespace utopia

#endif  // UTOPIA_KOKKOS_EVAL_REDUCE_HPP
