#ifndef UTOPIA_CUDA_BACKEND_TPL_HPP
#define UTOPIA_CUDA_BACKEND_TPL_HPP 

#include "utopia_CUDAMatrix.hpp"
#include "utopia_CUDAVector.hpp"

namespace utopia {
	namespace cuda_double {
		typedef utopia::CUDAMatrix<double> CUDAMatrix;
		typedef utopia::CUDAVector<double> CUDAVector;

		void describe(const CUDAMatrix &m);
		void describe(const CUDAVector &v);
		void build_identity(const int rows, const int cols, CUDAMatrix &m);	
		void build_values(const int n, const double value, CUDAVector &v);
                void build_values(const int n1, const int n2, const double value, CUDAMatrix &m);
		double dot(const CUDAVector &left, const CUDAVector &right);
		void mat_vec_mul(const CUDAMatrix &left, const CUDAVector &right, CUDAVector &result);
                void mat_mat_mul(const CUDAMatrix &left, const CUDAMatrix &right, const int width, CUDAMatrix &result);
	}	
}

#endif //UTOPIA_CUDA_BACKEND_TPL_HPP
