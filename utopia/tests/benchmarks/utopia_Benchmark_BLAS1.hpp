#ifndef UTOPIA_BENCHMARK_BLASS_1_HPP
#define UTOPIA_BENCHMARK_BLASS_1_HPP

#include "utopia_Chrono.hpp"
#include "utopia_MPI.hpp"
#include "utopia.hpp"
#include "utopia_Benchmark.hpp"
#include <string>

namespace utopia {
	
	template<class Matrix, class Vector>
	class BenchmarkBlas1 : public Benchmark {
	public:
		DEF_UTOPIA_SCALAR(Vector)

		virtual std::string name() override
		{
			return "blas1";
		}

		void initialize() override
		{
			const SizeType base_n = 500;
			const SizeType n_instances = 10;

			for(SizeType i = 0; i < n_instances; ++i) {
				const SizeType n = base_n * (i + 1);

				//axpy
				this->register_experiment(
					"axpy_" + std::to_string(i),
					[n]() {
						const Scalar alpha = 0.1;
						const Vector x = local_values(n, 1.);
						Vector y = local_values(n, 2.);
						
						y += alpha * x;
						y = alpha * x + y;
						y = y + alpha * x;
					}
				);


				//...
			}
		}

	};
}

#endif //UTOPIA_BENCHMARK_BLASS_1_HPP
