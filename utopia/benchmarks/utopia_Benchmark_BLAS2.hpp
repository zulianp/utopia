#ifndef UTOPIA_BENCHMARK_BLAS_2_HPP
#define UTOPIA_BENCHMARK_BLAS_2_HPP

#include "utopia_Chrono.hpp"
#include "utopia_MPI.hpp"
#include "utopia.hpp"
#include "utopia_Benchmark.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"

#include <string>

namespace utopia {
    //http://www.netlib.org/blas/#_level_2
    template<class Matrix, class Vector>
    class BenchmarkBlas2 : public Benchmark {
    public:
        using Traits   = utopia::Traits<Vector>;
        using Scalar   = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm     = typename Traits::Communicator;

        BenchmarkBlas2(const Comm &comm = Comm()) : comm_(comm) {}

        std::string name() override { return "BLAS 2"; }

        void initialize() override
        {
            static const bool is_sparse = utopia::is_sparse<Matrix>::value;

            const SizeType base_n = is_sparse? 9000 : 90;
            const SizeType n_instances = 10;

            for(SizeType i = 0; i < n_instances; ++i) {
                const SizeType n = base_n * (i + 1);

                auto vl = layout(comm_, n, n * comm_.size());
                auto ml = layout(comm_, n, n, n * comm_.size(), n * comm_.size());

                //measure allocation time of two vectors and the matrix
                this->register_experiment(
                    "allocation_" + std::to_string(i),
                    [vl, ml]() {
                        Vector x(vl, 1.);
                        Vector y(vl, 2.);
                        Matrix A; A.sparse(ml, 3, 2);
                    }
                );

                //measure assembly time of the operator
                this->register_experiment(
                    "assembly_" + std::to_string(i),
                    [vl, ml]() {
                        Matrix A; A.sparse(ml, 3, 2);
                        assemble_laplacian_1D(A);
                    }
                );

                //matrix-vector mult
                this->register_experiment(
                    "mv_" + std::to_string(i),
                    [vl, ml]() {
                        Vector x(vl, 1.);
                        Vector y(vl, 2.);
                        Matrix A; A.sparse(ml, 3, 2);
                        assemble_laplacian_1D(A);

                        y = A * x;
                    }
                );

                //...
            }
        }

        Comm comm_;

    };
}

#endif //UTOPIA_BENCHMARK_BLAS_2_HPP
