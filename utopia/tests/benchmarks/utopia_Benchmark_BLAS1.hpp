#ifndef UTOPIA_BENCHMARK_BLASS_1_HPP
#define UTOPIA_BENCHMARK_BLASS_1_HPP

#include "utopia_Chrono.hpp"
#include "utopia_MPI.hpp"
#include "utopia.hpp"
#include "utopia_Benchmark.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"
#include <string>
#include <cmath>
#include <cassert>

namespace utopia {
	//http://www.netlib.org/blas/#_level_1
	template<class Matrix, class Vector>
	class BenchmarkBlas1 : public Benchmark {
	public:
		DEF_UTOPIA_SCALAR(Vector)

		virtual std::string name() override
		{
			return "BLAS 1";
		}

		void initialize() override
		{
			const SizeType base_n = 9000;
			const SizeType n_instances = 10;

			for(SizeType i = 0; i < n_instances; ++i) {
				const SizeType n = base_n * (i + 1);
				//Vectors
					
				//measure allocation time of two vectors
				this->register_experiment(
					"vec_allocation_" + std::to_string(i),
					[n]() {
						Vector x = local_values(n, 1.);
						Vector y = local_values(n, 2.);
					}
				);
				
				//axpy
				this->register_experiment(
					"vec_axpy_" + std::to_string(i),
					[n]() {
						const Vector x = local_values(n, 1.);
						Vector y = local_values(n, 2.);
						const Scalar alpha = 0.1;
						y += alpha * x;
						y = alpha * x + y;
						y = y + alpha * x;
					}
				);

				//norms
				this->register_experiment(
					"vec_norms_" + std::to_string(i),
					[n]() {
						const Vector x = local_values(n, 1.);
						const Scalar norm2_x = norm2(x);
						const Scalar norm1_x = norm1(x);
						const Scalar norm_infty_x = norm_infty(x);

						//some testing
						assert(approxeq(norm2_x, Scalar(std::sqrt(n))));
						assert(approxeq(norm1_x, Scalar(n)));
						assert(approxeq(norm_infty_x, Scalar(1)));

					}
				);

				//scale
				this->register_experiment(
					"vec_scale_" + std::to_string(i),
					[n]() {
						Vector x = local_values(n, 1.);
						x *= 0.1;
						// x = x * 0.1;
					}
				);

				//dot
				this->register_experiment(
					"vec_dot_" + std::to_string(i),
					[n]() {
						const Vector x = local_values(n, 1.);
						const Vector y = local_values(n, 2.);
						const Scalar d = dot(x, y);

						assert(approxeq(d, Scalar(n*2)));
					}
				);

				//distance
				this->register_experiment(
					"vec_dist_" + std::to_string(i),
					[n]() {
						const Vector x = local_values(n, 1.);
						const Vector y = local_values(n, 2.);
						const Scalar d = norm2(x - y);

						assert(approxeq(d, Scalar(std::sqrt(n))));
					}
				);

				//Matrices
				//measure allocation time of one matrix
				this->register_experiment(
					"mat_allocation_" + std::to_string(i),
					[n]() {
						Matrix A = local_sparse(n, n, 3);
						assemble_laplacian_1D(A);
					}
				);

				//axpy
				this->register_experiment(
					"mat_axpy_" + std::to_string(i),
					[n]() {
						Matrix A = local_sparse(n, n, 3);
						assemble_laplacian_1D(A);

						Matrix B = A;
						Matrix C = A + B;
					}
				);

				//...
			}
		}

	};
}

#endif //UTOPIA_BENCHMARK_BLASS_1_HPP
