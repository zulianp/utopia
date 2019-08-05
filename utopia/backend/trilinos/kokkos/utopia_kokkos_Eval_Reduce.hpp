#ifndef UTOPIA_KOKKOS_EVAL_REDUCE_HPP
#define UTOPIA_KOKKOS_EVAL_REDUCE_HPP

#include "utopia_kokkos_Operations.hpp"
#include <Kokkos_Core.hpp>
#include <iostream>
#include <cassert>

namespace utopia {
    template<typename Data, typename KokkosOp, typename Scalar>
    struct OpFunctor {
        KOKKOS_INLINE_FUNCTION void join(volatile Scalar &val, const volatile Scalar &other) const {
            // Kokkos forces us to have the input values being declared volatile. Hence we need to make copies for the reduction operations
            double tmp1=val, tmp2=other;
            val = _op.apply(tmp1, tmp2);
        }
        KOKKOS_INLINE_FUNCTION void operator()(const int& i, Scalar &val) const {
            val = _op.apply(val, _data(i, 0));
        }
        const KokkosOp& _op;
        const Data& _data;
    };


    template<class Vector, class Op>
    class KokkosEvalReduce {
    public:
        inline static typename Vector::Scalar eval(const Vector &vec, const Op op, const typename Vector::Scalar &initial_value)
        {
            using ExecutionSpaceT = typename Vector::vector_type::execution_space;
            using Scalar = typename Vector::Scalar;
            using Data = decltype(vec.implementation().template getLocalView<ExecutionSpaceT>());

            assert(!vec.empty());
            auto data = vec.implementation().template getLocalView<ExecutionSpaceT>();

            Scalar ret = initial_value;
            KokkosOp<Scalar, Op> k_op;
            OpFunctor<Data, KokkosOp<Scalar, Op>, Scalar> functor{k_op, data};
            Kokkos::parallel_reduce(data.extent(0), functor, ret);

            return ret;
        }
    };

}

#endif //UTOPIA_KOKKOS_EVAL_REDUCE_HPP
