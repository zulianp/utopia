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

    template<class Matrix, class Vector>
    class BenchmarkTransform : public Benchmark {
    public:
        using Traits   = utopia::Traits<Vector>;
        using Scalar   = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm     = typename Traits::Communicator;

        BenchmarkTransform(const Comm &comm = Comm()) : comm_(comm) {}


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

                auto vl = layout(comm_, 5 * n, 5 * n * comm_.size());
                auto ml = layout(comm_, n, n, n * comm_.size(), n * comm_.size());

                this->register_experiment(
                    "vec_transform_" + std::to_string(i),
                    [vl]() {
                        Vector v(vl, -3.0);

                        UTOPIA_NO_ALLOC_BEGIN("vec_transform");

                        v = abs(v);
                        v = sqrt(v);
                        v = pow2(v);
                        v = -v;
                        v = exp(v);

                        UTOPIA_NO_ALLOC_END();
                    }
                );

                this->register_experiment(
                    "mat_transform_" + std::to_string(i),
                    [ml]() {
                        Matrix A; A.sparse(ml, 3, 2);
                        assemble_laplacian_1D(A);

                        // UTOPIA_NO_ALLOC_BEGIN("mat_transform");

                        A = abs(A);
                        A = sqrt(A);
                        A = pow2(A);
                        A = -A;
                        A = exp(A);

                        // UTOPIA_NO_ALLOC_END();
                    }
                );

                this->register_experiment(
                    "mat_chop_" + std::to_string(i),
                    [ml]() {
                        Matrix A; A.sparse(ml, 3, 2);
                        assemble_laplacian_1D(A);

                        chop(A, 0.001);
                        chop_smaller_than(A, 0.001);
                        chop_greater_than(A, 2.0);
                    }
                );
            }
        }

    private:
      Comm comm_;

    };
}

#endif //UTOPIA_BENCHMARK_TRANSFORM_HPP
