#ifndef UTOPIA_BENCHMARK_ACCESS_HPP
#define UTOPIA_BENCHMARK_ACCESS_HPP

#include <cassert>
#include <cmath>
#include <string>
#include "utopia_assemble_laplacian_1D.hpp"
#include "utopia.hpp"
#include "utopia_Benchmark.hpp"
#include "utopia_Chrono.hpp"
#include "utopia_MPI.hpp"
#include "utopia_RangeDevice.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class BenchmarkAccess : public Benchmark {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        BenchmarkAccess(const Comm &comm = Comm()) : comm_(comm) {}

        std::string name() override { return "Access"; }

        void initialize() override {
            static const bool is_sparse = utopia::is_sparse<Matrix>::value;

            const SizeType base_n = is_sparse ? 9000 : 90;
            const SizeType n_instances = 10;

            for (SizeType i = 0; i < n_instances; ++i) {
                const SizeType n = base_n * (i + 1);

                auto vl = layout(comm_, n, n * comm_.size());
                auto ml = layout(comm_, n, n, n * comm_.size(), n * comm_.size());

                // Vectors
                this->register_experiment("vec_set_" + std::to_string(i), [vl]() {
                    Vector x(vl, 1.);
                    auto r = range(x);

                    SizeType K = 100;
                    for (SizeType k = 0; k < K; ++k) {
                        Write<Vector> w_(x);
                        for (auto i = r.begin(); i < r.end(); ++i) {
                            x.set(i, 2.);
                        }
                    }

                    utopia_test_assert(approxeq(sum(x), size(x).get(0) * 2.));
                });

                this->register_experiment("vec_get_" + std::to_string(i), [vl]() {
                    Vector x(vl, 1.);

                    auto r = range(x);

                    Scalar res = 0.;
                    SizeType K = 100;

                    for (SizeType k = 0; k < K; ++k) {
                        Read<Vector> w_(x);
                        for (auto i = r.begin(); i < r.end(); ++i) {
                            res += x.get(i) / K;
                        }
                    }

                    utopia_test_assert(approxeq(sum(x), size(x).get(0) * 1.));
                });

                this->register_experiment("vec_set_get_" + std::to_string(i), [vl]() {
                    Vector x(vl, 1.);
                    Vector y(vl, 0.0);

                    auto r = range(x);

                    {
                        Read<Vector> r_(x);
                        Write<Vector> w_(y);

                        for (auto i = r.begin(); i < r.end(); ++i) {
                            y.set(i, x.get(i));
                        }
                    }

                    utopia_test_assert(approxeq(sum(y), size(x).get(0) * 1.));
                });

                // measure loop time for vectors
                this->register_experiment("vec_for_" + std::to_string(i), [vl]() {
                    Vector x(vl, 1.);

                    {
                        auto x_view = view_device(x);
                        parallel_for(
                            range_device(x), UTOPIA_LAMBDA(const SizeType i) { x_view.set(i, i); });
                    }

                    Scalar res = sum(x);

                    res /= size(x).get(0);

                    {
                        auto x_view = local_view_device(x);
                        parallel_for(
                            local_range_device(x),
                            UTOPIA_LAMBDA(const SizeType i) { x_view.set(i, x_view.get(i) - res); });
                    }
                });

                // Matrices
                this->register_experiment("mat_assemble_lapl_" + std::to_string(i), [ml]() {
                    Matrix A;
                    A.sparse(ml, 3, 2);
                    assemble_laplacian_1D(A);
                });

                this->register_experiment("mat_each_read_" + std::to_string(i), [ml]() {
                    Matrix A;
                    A.sparse(ml, 3, 2);
                    assemble_laplacian_1D(A);
                    A.read(UTOPIA_LAMBDA(const SizeType /*i*/, const SizeType /*j*/, const Scalar){});
                });

                //...
            }
        }

    private:
        Comm comm_;
    };
}  // namespace utopia

#endif  // UTOPIA_BENCHMARK_ACCESS_HPP
