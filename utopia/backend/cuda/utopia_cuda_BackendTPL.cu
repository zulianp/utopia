#include "utopia_cuda_BackendTPL.hpp"

#include <thrust/inner_product.h>
#include <utility>

#include "utopia_cuda_Error.hpp"

namespace utopia {

	namespace cuda_generic {
		template<typename T>
		__global__ void build_identity(const int rows, const int cols, T *values)
		{
			int id_x = blockIdx.x * blockDim.x + threadIdx.x;
			int id_y = blockIdx.y * blockDim.y + threadIdx.y;

			if(id_x < rows && id_y < cols) {
				if(id_x == id_y) {
					values[id_x * cols + id_y] = 1;
				} else {
					values[id_x * cols + id_y] = 0;
				}
			}
		}

		template<typename T>
		__global__ void mat_vec_mul(const int rows, const int cols, const T *mat_left, const T *vec_right, T *result)
		{
			int id_x = blockIdx.x * blockDim.x + threadIdx.x;
			if(id_x < rows) {
				
				const int offset_i = id_x * cols;

				T prod = 0;
				for(int j = 0; j < cols; ++j) {
					prod += mat_left[offset_i + j] * vec_right[j];
				}

				// __synchthreads();
				result[id_x] = prod;
			}
		}

		std::pair<dim3, dim3> get_sizes_2(const int n_x, const int n_y)
		{
			// dim3 n_blocks( ceil(n_x/22.), ceil(n_y/22.) );
			// dim3 n_threads_x_block(n_x/n_blocks.x + 1, 
			// 					   n_y/n_blocks.y + 1);

			dim3 n_blocks(n_x, n_y);
			dim3 n_threads_x_block(1, 1);

			return std::make_pair(n_blocks, n_threads_x_block);
		}

		std::pair<dim3, dim3> get_sizes_1(const int n_x)
		{
			dim3 n_blocks(ceil(n_x/512.));
			dim3 n_threads_x_block(ceil(double(n_x)/n_blocks.x));
			return std::make_pair(n_blocks, n_threads_x_block);
		}
	}


	namespace cuda_double {

		template<class Tensor>
		static void describe_values(const Tensor &t) 
		{
			thrust::copy(t.values.begin(), t.values.end(), std::ostream_iterator<double>(std::cout, "\n"));
		}

		void describe(const CUDAMatrix &m)
		{
			describe_values(m);
		}

		void describe(const CUDAVector &v)
		{
			describe_values(v);
		}

		void build_identity(const int rows, const int cols, CUDAMatrix &m)
		{
			m.values.resize(rows * cols);
			m.rows = rows;
			m.cols = cols;

			double *raw_ptr = thrust::raw_pointer_cast(&m.values[0]);
			dim3 n_blocks, n_threads_x_block;
			std::pair<dim3, dim3> s = cuda_generic::get_sizes_2(rows, cols);

			cuda_generic::build_identity<double><<<s.first, s.second>>>(rows, cols, raw_ptr);
		}
		
		void build_values(const int n, const double value, CUDAVector &v)
		{
			v.values.resize(n);
			thrust::fill(v.values.begin(), v.values.end(), value);
		}

		double dot(const CUDAVector &left, const CUDAVector &right)
		{
			return thrust::inner_product(left.values.begin(), left.values.end(), right.values.begin(), 0.0);
		}

		void mat_vec_mul(const CUDAMatrix &left, const CUDAVector &right, CUDAVector &result)
		{
			result.values.resize(left.rows);

			dim3 n_blocks, n_threads_x_block;
			
			std::pair<dim3, dim3> s = cuda_generic::get_sizes_1(left.rows);

			cuda_generic::mat_vec_mul<double><<<s.first, s.second>>>(
				left.rows,
				left.cols,
				thrust::raw_pointer_cast(&left.values[0]),
				thrust::raw_pointer_cast(&right.values[0]),
				thrust::raw_pointer_cast(&result.values[0])
				);

			CUDAError::CheckLastError();
		}
	}
}

