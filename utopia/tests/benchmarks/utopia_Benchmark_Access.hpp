#ifndef UTOPIA_BENCHMARK_ACCESS_HPP
#define UTOPIA_BENCHMARK_ACCESS_HPP

#include "utopia_Chrono.hpp"
#include "utopia_MPI.hpp"
#include "utopia.hpp"
#include "utopia_Benchmark.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"
#include <string>
#include <cmath>
#include <cassert>

namespace utopia {

	template<class Matrix, class Vector>
	class BenchmarkAccess : public Benchmark {
	public:
		DEF_UTOPIA_SCALAR(Vector)
		using SizeType = UTOPIA_SIZE_TYPE(Vector);

		virtual std::string name() override
		{
			return "Access";
		}

		void initialize() override
		{
			const SizeType base_n = 9000;
			const SizeType n_instances = 10;

			for(SizeType i = 0; i < n_instances; ++i) {
				const SizeType n = base_n * (i + 1);
				//Vectors
					
				//measure loop time for vectors
				this->register_experiment(
					"vec_each" + std::to_string(i),
					[n]() {
						Vector x = local_values(n, 1.);

						each_write(x, [](const SizeType i) -> Scalar {
							return i;
						});

						Scalar res = 0.0;
						each_read(x, [&res](const SizeType i, const Scalar val) {
							res += val;
						});

						res /= size(x).get(0);

						each_transform(x, x, [res](const SizeType i, const Scalar val) -> Scalar {
							return val - res;
						});
					}
				);

				this->register_experiment(
					"mat_each_read" + std::to_string(i),
					[n]() {
						
						Matrix A = local_sparse(n, n, 3);
						assemble_laplacian_1D(A);

						auto N = size(A).get(0);

						Scalar res = 0.0;
						each_read(A, [&res](const SizeType i, const SizeType j, const Scalar val) {
							res += val;
						});

						utopia_test_assert(approxeq(res, 0.));
					}
				);
				
				//...
			}
		}

	};
}

#endif //UTOPIA_BENCHMARK_ACCESS_HPP
