#ifndef UTOPIA_CUDA_VECTOR_HPP
#define UTOPIA_CUDA_VECTOR_HPP 

#include <thrust/device_vector.h>

namespace utopia {
	template<typename _Scalar>
	class CUDAVector {
	public:
		typedef _Scalar Scalar;
		thrust::device_vector<Scalar> values;
	};
}

#endif