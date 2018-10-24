#ifndef UTOPIA_BENCHMARK_ALGORITHMS_HPP
#define UTOPIA_BENCHMARK_ALGORITHMS_HPP

#include "utopia_Chrono.hpp"
#include "utopia_MPI.hpp"
#include "utopia.hpp"
#include "utopia_Benchmark.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"
#include "test_problems/utopia_RastriginTestFunction.hpp"
#include "utopia_Preconditioner.hpp"
#include "utopia_ProjectedGradient.hpp"
#include "utopia_ProjectedConjugateGradient.hpp"
#include "utopia_ProjectedGaussSeidel.hpp"

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
						auto A_op = utopia::op_ref(A);
						cg.solve(*A_op, b, x);

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


				this->register_experiment(
					"projected_gradient" + std::to_string(i),
					[i]() {
						ProjectedGradient<Matrix, Vector, HOMEMADE> pg;
						run_qp_solver((base_n/2) * (i + 1), pg);
					}
				);

				this->register_experiment(
					"projected_conjugate_gradient" + std::to_string(i),
					[i]() {
						ProjectedConjugateGradient<Matrix, Vector, HOMEMADE> pg;
						run_qp_solver((base_n/2) * (i + 1), pg);
					}
				);

				this->register_experiment(
					"projected_gauss_seidel" + std::to_string(i),
					[i]() {
						ProjectedGaussSeidel<Matrix, Vector, HOMEMADE> pg;
						run_qp_solver((base_n/2) * (i + 1), pg);
					}
				);

			}
		}

	private:

		template<class QPSolver>
		static void run_qp_solver(const SizeType n, QPSolver &qp_solver) {

		    Matrix m = local_sparse(n, n, 3);
		    assemble_laplacian_1D(m);

		    auto N = size(m).get(0);

		    {
		        Range r = row_range(m);
		        Write<Matrix> w(m);
		        if(r.begin() == 0) {
		            m.set(0, 0, 1.);
		            m.set(0, 1, 0);
		        }

		        if(r.end() == N) {
		            m.set(N-1, N-1, 1.);
		            m.set(N-1, N-2, 0);
		        }
		    }

		    Vector rhs = local_values(n, 1.);
		    {
		        //Creating test vector (alternative way see [assemble vector alternative], which might be easier for beginners)
		        Range r = range(rhs);
		        Write<Vector> w(rhs);

		        if(r.begin() == 0) {
		            rhs.set(0, 0);
		        }

		        if(r.end() == N) {
		            rhs.set(N-1, 0.);
		        }
		    }

		    Vector upper_bound = local_values(n, 100.0);
		    Vector solution    = local_zeros(n);


		    qp_solver.max_it(N*2);
		    // qp_solver.verbose(true);
		    qp_solver.set_box_constraints(make_upper_bound_constraints(make_ref(upper_bound)));

		    bool ok = qp_solver.solve(m, rhs, solution);
		    assert(ok);
		}

	};
}

#endif //UTOPIA_BENCHMARK_ALGORITHMS_HPP
