#ifndef UTOPIA_BENCHMARK_TRANSFORM_HPP
#define UTOPIA_BENCHMARK_TRANSFORM_HPP

#include "utopia_Chrono.hpp"
#include "utopia_MPI.hpp"
#include "utopia.hpp"
#include "utopia_Benchmark.hpp"
#include "utopia_MatChop.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"
#include <string>
#include <cmath>
#include <cassert>

namespace utopia {
    //http://www.netlib.org/blas/#_level_1
    template<class Matrix, class Vector>
    class BenchmarkTransform : public Benchmark {
    public:
        DEF_UTOPIA_SCALAR(Vector);

        virtual std::string name() override
        {
            return "Transform";
        }

        void initialize() override
        {
            static const bool is_sparse = utopia::is_sparse<Matrix>::value;

            const SizeType base_n = is_sparse? 9000 : 90;
            const SizeType n_instances = 10;

            for(SizeType i = 0; i < n_instances; ++i) {
                const SizeType n = base_n * (i + 1);
              
                this->register_experiment(
                    "mat_transform_" + std::to_string(i),
                    [n]() {
                        Matrix A = local_sparse(n, n, 3);
                        assemble_laplacian_1D(A);

                        A = abs(A);
                        A = sqrt(A);
                        A = pow2(A);
                        A = -A;
                        A = exp(A);
                    }
                );

                this->register_experiment(
                    "mat_chop_" + std::to_string(i),
                    [n]() {
                        Matrix A = local_sparse(n, n, 3);
                        assemble_laplacian_1D(A);

                        chop(A, 0.001);
                        chop_smaller_than(A, 0.001);
                        chop_greater_than(A, 2.0);
                    }
                );
            }
        }

    };
}

#endif //UTOPIA_BENCHMARK_TRANSFORM_HPP