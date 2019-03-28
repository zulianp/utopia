#ifndef UTOPIA_BENCHMARK_BLAS_3_HPP
#define UTOPIA_BENCHMARK_BLAS_3_HPP

#include "utopia_Chrono.hpp"
#include "utopia_MPI.hpp"
#include "utopia.hpp"
#include "utopia_Benchmark.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"

#include <string>

namespace utopia {
    //http://www.netlib.org/blas/#_level_3
    template<class Matrix, class Vector>
    class BenchmarkBlas3 : public Benchmark {
    public:
        DEF_UTOPIA_SCALAR(Vector)

        virtual std::string name() override
        {
            return "BLAS 3";
        }

        void initialize() override
        {
            const SizeType base_n = 9000;
            const SizeType n_instances = 10;

            for(SizeType i = 0; i < n_instances; ++i) {
                const SizeType n = base_n * (i + 1);

                //measure allocation time of two vectors and the matrix
                this->register_experiment(
                    "allocation_" + std::to_string(i),
                    [n]() {
                        Matrix A = local_sparse(n, n, 3);
                    }
                );

                //measure assembly time of the operator
                this->register_experiment(
                    "assembly_" + std::to_string(i),
                    [n]() {
                        Matrix A = local_sparse(n, n, 3);
                        assemble_laplacian_1D(A);
                        //matrix copy
                        Matrix B = A;
                    }
                );

                //matrix-vector mult
                this->register_experiment(
                    "mat_mat_mult_" + std::to_string(i),
                    [n]() {
                        Matrix A = local_sparse(n, n, 3);
                        assemble_laplacian_1D(A);
                        //matrix copy
                        Matrix B = A;
                        Matrix C = A * B;
                    }
                );

                //...
            }
        }

    };
}

#endif //UTOPIA_BENCHMARK_BLAS_3_HPP
