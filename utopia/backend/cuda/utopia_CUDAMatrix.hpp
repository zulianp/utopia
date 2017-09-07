#ifndef UTOPIA_CUDA_MATRIX_HPP
#define UTOPIA_CUDA_MATRIX_HPP 

#include <thrust/device_vector.h>

namespace utopia {
	template<typename _Scalar>
	class CUDAMatrix {
	public:
		typedef _Scalar Scalar;

		thrust::device_vector<Scalar> values;
		int rows, cols;
	};
}

#endif
