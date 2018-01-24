#ifndef UTOPIA_CUDA_BACKEND_TPL_HPP
#define UTOPIA_CUDA_BACKEND_TPL_HPP 

#include "utopia_cuda_Matrix.hpp"
#include "utopia_cuda_Vector.hpp"

namespace utopia {
	namespace cuda_double {
		typedef utopia::CUDAMatrix<double> CUDAMatrix;
		typedef utopia::CUDAVector<double> CUDAVector;

		void describe(const CUDAMatrix &m);
		void describe(const CUDAVector &v);
		void build_identity(const int rows, const int cols, CUDAMatrix &m);	
		void build_values(const int n, const double value, CUDAVector &v);
		double dot(const CUDAVector &left, const CUDAVector &right);
		void mat_vec_mul(const CUDAMatrix &left, const CUDAVector &right, CUDAVector &result);
	}	
}

#endif //UTOPIA_CUDA_BACKEND_TPL_HPP
