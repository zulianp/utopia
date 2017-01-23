#ifndef UTOPIA_CUDA_TRAITS_HPP
#define UTOPIA_CUDA_TRAITS_HPP

#include "utopia_Traits.hpp"
#include "utopia_CUDAMatrix.hpp"
#include "utopia_CUDAVector.hpp"

namespace utopia {
    template<typename T>
	class CudaTraits {
	public:
		typedef T Scalar;
		typedef utopia::CUDAMatrix<T> Matrix;
		typedef utopia::CUDAVector<T> Vector;
		typedef long SizeType;
		enum {
			Backend = CUDA
		};
	};

	UTOPIA_MAKE_TRAITS_TPL_1(CUDAVector, CudaTraits);
	UTOPIA_MAKE_TRAITS_DENSE_TPL_1(CUDAMatrix, CudaTraits);
}

#endif //UTOPIA_CUDA_TRAITS_HPP


