#ifndef UTOPIA_BENCHMARK_ALGORITHMS_HPP
#define UTOPIA_BENCHMARK_ALGORITHMS_HPP

#include "utopia_Chrono.hpp"
#include "utopia_MPI.hpp"
#include "utopia.hpp"
#include "utopia_Benchmark.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"
#include "test_problems/utopia_RastriginTestFunction.hpp"

#include <string>
#include <cassert>

namespace utopia {

	template<class Matrix, class Vector>
	class BenchmarkAlgorithms : public Benchmark {
	public:
		DEF_UTOPIA_SCALAR(Vector)

		virtual std::string name() override
		{
			return "Algorithms";
		}

		void initialize() override
		{
			const SizeType base_n = 1000;
			const SizeType n_instances = 5;

			for(SizeType i = 0; i < n_instances; ++i) {
				const SizeType n = base_n * (i + 1);
					
				//Conjugate gradient method
				this->register_experiment(
					"cg_" + std::to_string(i),
					[n]() {
						Matrix A = local_sparse(n, n, 3); 
						Vector b = local_values(n, 1.);
						Vector x = local_values(n, 0.);

						assemble_laplacian_1D(A, true);

						auto N = size(A).get(0);
						{
							Range r = row_range(A);
							Write<Vector> w_b(b);
							
							if(r.begin() == 0) {
								b.set(0, 0.);
							}

							if(r.end() == N) {
								b.set(N-1, 0.);
							}
						}

						ConjugateGradient<Matrix, Vector, HOMEMADE> cg;
						cg.max_it(N);
						cg.solve(A, b, x);

						double err = norm2(b - A * x);
						assert(approxeq(A * x, b, 1e-6));
					}
				);	

				this->register_experiment(
					"newton_cg_" + std::to_string(i),
					[i]() {
						Rastrigin<Matrix, Vector> fun;
						Vector x = local_values(10 * (i+1), 1.);

						ConjugateGradient<Matrix, Vector, HOMEMADE> cg;
						cg.max_it(size(x).get(0));

						auto backtracking = std::make_shared<utopia::Backtracking<Matrix, Vector> >();

						Newton<Matrix, Vector, HOMEMADE> newton(make_ref(cg));
						newton.set_line_search_strategy(backtracking);

						double mag_x0 = -1;
						fun.value(x, mag_x0);

						newton.solve(fun, x);

						double mag_x = -1.;
						fun.value(x, mag_x);
						assert(mag_x <= mag_x0);
					}
				);

				this->register_experiment(
					"trust_region_" + std::to_string(i),
					[i]() {
						Rastrigin<Matrix, Vector> fun;
						Vector x = local_values(10 * (i+1), 1.);

						TrustRegion<Matrix, Vector> trust_region;
						trust_region.verbose(false); 

						double mag_x0 = -1;
						fun.value(x, mag_x0);

						trust_region.solve(fun, x);

						double mag_x = -1.;
						fun.value(x, mag_x);
						assert(mag_x <= mag_x0);
					}
				);

			}
		}

	};
}

#endif //UTOPIA_BENCHMARK_ALGORITHMS_HPP
