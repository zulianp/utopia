#ifndef UTOPIA_BENCHMARK_BLASS_1_HPP
#define UTOPIA_BENCHMARK_BLASS_1_HPP

#include "utopia_Chrono.hpp"
#include "utopia_MPI.hpp"
#include "utopia.hpp"
#include "utopia_Benchmark.hpp"
#include <string>

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
					
				//measure allocation time of two vectors
				this->register_experiment(
					"allocation_" + std::to_string(i),
					[n]() {
						Vector x = local_values(n, 1.);
						Vector y = local_values(n, 2.);
					}
				);
				
				//axpy
				this->register_experiment(
					"axpy_" + std::to_string(i),
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
					"norms_" + std::to_string(i),
					[n]() {
						const Vector x = local_values(n, 1.);
						const Scalar norm2_x = norm2(x);
						const Scalar norm1_x = norm1(x);
						const Scalar norm_infty_x = norm_infty(x);
					}
				);


				//scale
				this->register_experiment(
					"scale_" + std::to_string(i),
					[n]() {
						Vector x = local_values(n, 1.);
						x *= 0.1;
						// x = x * 0.1;
					}
				);


				//dot
				this->register_experiment(
					"dot_" + std::to_string(i),
					[n]() {
						const Vector x = local_values(n, 1.);
						const Vector y = local_values(n, 2.);
						const Scalar d = dot(x, y);
					}
				);
				

				//...
			}
		}

	};
}

#endif //UTOPIA_BENCHMARK_BLASS_1_HPP
