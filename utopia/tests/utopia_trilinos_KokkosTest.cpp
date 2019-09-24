
#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS

#include "utopia_trilinos_KokkosTest.hpp"
#include "utopia.hpp"
#include "utopia_trilinos.hpp"
#include "utopia_trilinos_solvers.hpp"
#include "utopia_kokkos_ParallelEach.hpp"
#include "utopia_trilinos_Each_impl.hpp"
#include <algorithm>

#include "utopia_Structure.hpp"
#include "utopia_Eval_Structure.hpp"
#include "test_problems/utopia_BratuMultilevelTestProblem.hpp"
#include "test_problems/utopia_TestProblems.hpp"

#include "utopia_IPTransfer.hpp"
#include <cmath>

namespace utopia {

    void kokkos_max()
    {
        auto n = 10;

        TpetraVectord v = local_values(n, 1.);

        TpetraVectord w = local_values(n, 10.);

        TpetraVectord z = max(w, v);

        utopia_test_assert(approxeq(z, w));
    }

    void kokkos_min()
    {
        auto n = 10;

        TpetraVectord v = local_values(n, 1.);

        TpetraVectord w = local_values(n, 10.);

        TpetraVectord z = min(w , v);

        utopia_test_assert(approxeq(z, v));
    }

    void kokkos_sum_reduction()
    {
        auto n = 10;

        TpetraVectord v = local_values(n, 1.);

        double z = sum(v);

        utopia_test_assert(approxeq(z, size(v).get(0) * 1.));
    }

    void kokkos_min_reduction()
    {
        auto n = 10;

        TpetraVectord v = local_values(n, 1.);

        auto r = range(v);

        each_write(v,  KOKKOS_LAMBDA(const SizeType i) -> double {
            return i - r.begin();
        });

        double z = min(v);

        utopia_test_assert(approxeq(0., z));
    }

    void kokkos_max_reduction()
    {
        auto n = 10;

        TpetraVectord v = local_values(n, 1.);

        auto r = range(v);

        each_write(v, KOKKOS_LAMBDA(const SizeType i) -> double {
            return i - r.begin();
        });

        double z = max(v);

        utopia_test_assert(approxeq(9., z));
    }

    void kokkos_write() {
        auto n = 10;

        TpetraVectord w = local_values(n, -1.);

        parallel_each_write(w, KOKKOS_LAMBDA(const SizeType i) -> double {
            return i;
        });

        {
            Read<TpetraVectord> r_(w);
            auto r = range(w);

            for(auto i = r.begin(); i < r.end(); ++i) {
                utopia_test_assert(approxeq(w.get(i), double(i)));
            }
        }
    }

    void kokkos_parallel_each_mat() {
        auto n = 10;

        TpetraMatrixd w = local_identity(n, n);

        parallel_each_write(w, KOKKOS_LAMBDA(const SizeType i, const SizeType j) -> double {
            return i * n + j;
        });

        //serial implementation for test
        each_read(w, [=](const SizeType i, const SizeType j, const double val) {
            utopia_test_assert(approxeq(double(i * n + j), val));
        });
    }

    void kokkos_read() {

        auto n = 10;
        TpetraVectord w = local_values(n, 50);

        parallel_each_read(w, KOKKOS_LAMBDA(const SizeType i, const double entry)
        { });
    }

    void kokkos_apply()
    {
        auto nr = 3;
        auto nc = 3;

        TpetraMatrixd P = local_sparse(nr, nc, 1);

        {
            Write<TpetraMatrixd> w_(P);
            auto r = row_range(P);
            auto cols = size(P).get(1);
            for(auto i = r.begin(); i < r.end(); ++i) {
                if(i >= cols) {
                    break;
                }

                P.set(i, i, 1.);

            }
        }

        parallel_each_apply(P, KOKKOS_LAMBDA(const double value) -> double {
            return value * 2.;
        });

    }

    void run_kokkos_test()
    {
        UTOPIA_UNIT_TEST_BEGIN("KokkosTest");
        UTOPIA_RUN_TEST(kokkos_max);
        UTOPIA_RUN_TEST(kokkos_min);
        UTOPIA_RUN_TEST(kokkos_min_reduction);
        UTOPIA_RUN_TEST(kokkos_max_reduction);
        UTOPIA_RUN_TEST(kokkos_sum_reduction);
        UTOPIA_RUN_TEST(kokkos_write);
        UTOPIA_RUN_TEST(kokkos_read);
        UTOPIA_RUN_TEST(kokkos_apply);
        UTOPIA_RUN_TEST(kokkos_parallel_each_mat);
        UTOPIA_UNIT_TEST_END("KokkosTest");
    }
}

#else //WITH_TRILINOS

namespace utopia
{
    void run_kokkos_test() {}
}

#endif //WITH_TRILINOS