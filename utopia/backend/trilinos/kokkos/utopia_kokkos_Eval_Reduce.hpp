#ifndef UTOPIA_KOKKOS_EVAL_REDUCE_HPP
#define UTOPIA_KOKKOS_EVAL_REDUCE_HPP

#include "utopia_kokkos_Operations.hpp"
#include <Kokkos_Core.hpp>
#include <iostream>
#include <cassert>

namespace utopia {

	template<class Vector, class Op> 
	class KokkosEvalReduce {
	public:
		inline static typename Vector::Scalar eval(const Vector &vec, const Op op, const typename Vector::Scalar &initial_value)
		{
			using ExecutionSpaceT = typename Vector::vector_type::execution_space;
			using Scalar = typename Vector::Scalar;
		    
		    assert(!vec.empty());
		    auto data = vec.implementation().template getLocalView<ExecutionSpaceT>();

            Scalar ret = initial_value;

            KokkosOp<Scalar, Op> k_op;
            
            Kokkos::parallel_reduce(data.extent(0), KOKKOS_LAMBDA(const int i, Scalar &val) {
                val = k_op.apply(val, data(i, 0));
            }, ret);

    	    return ret;
		}
	};

}

#endif //UTOPIA_KOKKOS_EVAL_REDUCE_HPP
