//
// Created by Eric Botter on 11/08/17.
//

#ifndef UTOPIA_UTOPIA_BENCHMARK_H
#define UTOPIA_UTOPIA_BENCHMARK_H

#include "test_problems/utopia_TestProblems.hpp"
#include "utopia.hpp"

namespace utopia {

	template <typename Matrix, typename SparseMatrix, typename Vector>
	class Benchmark {
	public:

		explicit Benchmark(size_t size) {
			size_ = size;
			sink_ = 0;
		}

		// call this to forbid compilers to optimize away computations
		void markUsed(void* value) {
			sink_ = *((char*)value);
		}

		void vectorAlgebra() {
			Vector v = zeros(size_);
			{
				Write<Vector> w(v);
				for (int i = 0; i < size_; ++i) {
					v.set(i, 1.5 * i);
				}
			}

			double n = norm2(v);
			n += norm_infty(v);
			v *= 2.5;
			n += norm2(v);
			n += norm_infty(v);
			n += dot(v, v * 0.1);
			n += norm2(v * (1.0 / double(norm2(v))));
			n += norm2(v) * norm2(v) / dot(v, v);

			markUsed(&n);
		}

		void matrixAlgebra() {
			size_t dim = static_cast<size_t>(std::sqrt(size_));
			Matrix m = identity(dim, dim);
			{
				Write<Matrix> w(m);
				for (int i = 0; i < dim; ++i) {
					for (int j = 0; j < dim; ++j) {
						m.set(i, j, i * 1.5 + j * 2.5);
					}
				}
			}
			Vector v1 = values(dim, 0.85);

			Vector v = m * v1;
			v1 = transpose(m) * v;
			v = m * v1;
			double n = sum(v);

			markUsed(&n);
		}

		void sparseMatrixAlgebra() {
			size_t dim = static_cast<size_t>(std::sqrt(size_));

			Vector v = values(dim, 1.2);
			SparseMatrix s1 = diag(v);
			SparseMatrix s2 = identity(dim, dim);
			{
				Write<SparseMatrix> w(s2);
				s2.set(0, 1, 1);
			}

			SparseMatrix s3 = s1 * s2 * transpose(s1) * s2 * transpose(s2);
			s1 = s3 + s2 - s1;
			s3 = s1 * s2 + s3 * transpose(s2);

			markUsed(&s3);
		}

		void mixedAlgebra() {
			size_t dim = static_cast<size_t>(std::sqrt(size_));
			Matrix m  = values(dim, dim, 0.001);
			Vector v1 = values(dim, 0.1);
			Vector v2 = values(dim, 0.2);
			Vector v3 = values(dim, 0.3);
			SparseMatrix s1 = diag(v1);

			Vector v = abs(m * v1 + v2 - v3);
			double n = fabs(sum(v));
			v = 10. * v1 + 5. * v2;
			n += fabs(sum(v));
			v = (0.1 * m - m) * v1;
			n += sum(v);
			v = abs(m * pow2(v1) + sqrt(v2) - v3);
			n += fabs(sum(v));
			v = abs(m * (transpose(0.1 * m - m) * (v1)));
			n += sum(v);
			v = (0.1 * pow2(v1) - pow2(v2)) + abs(v3);
			n += fabs(sum(v));
			v = m * v2 * 0.1;
			n += fabs(sum(v));
			v = (m + 0.1 * identity(dim, dim)) * values(dim, 0.5);
			n += fabs(sum(v));
			v = s1 * v - m * v2 + s1 * (v3 - v1);
			n += fabs(sum(v));

			markUsed(&n);
		}

		void solvers() {
			size_t dim = static_cast<size_t>(std::sqrt(size_));

			Vector v = values(dim, 1.2);
			TestFunctionND_1<Matrix, Vector> fun(static_cast<int>(dim));
			solve(fun, v);

			auto linear_solver  = std::make_shared< ConjugateGradient<Matrix, Vector> >();
			auto preconditioner = std::make_shared< InvDiagPreconditioner<Matrix, Vector> >();
			linear_solver->set_preconditioner(preconditioner);
			Newton<Matrix, Vector> newton_solver(linear_solver);
			v += values(dim, 0.8);
			newton_solver.solve(fun, v);

			markUsed(&v);
		}

		void runAll() {
			vectorAlgebra();
			matrixAlgebra();
			if (mpi_world_size() == 1) {
				sparseMatrixAlgebra();
			}
			mixedAlgebra();
			solvers();
		}

	private:
		size_t size_;
		volatile char sink_;
	};
}


#endif //UTOPIA_UTOPIA_BENCHMARK_H
