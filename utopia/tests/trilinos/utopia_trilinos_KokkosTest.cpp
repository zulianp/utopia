
#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS
#include <algorithm>
#include "utopia.hpp"
#include "utopia_Testing.hpp"
//#include "utopia_kokkos_ParallelEach.hpp"
#include "utopia_trilinos.hpp"
#include "utopia_trilinos_Each_impl.hpp"
#include "utopia_trilinos_solvers.hpp"

#include "test_problems/utopia_TestProblems.hpp"
#include "utopia_Eval_Structure.hpp"
#include "utopia_Structure.hpp"

#include <cmath>
#include "utopia_IPTransfer.hpp"

namespace utopia {

    class KokkosTest {
    public:
        using Traits = utopia::Traits<TpetraVectord>;
        using Comm = Traits::Communicator;
        using SizeType = Traits::SizeType;
        using Scalar = Traits::Scalar;

        void kokkos_max() {
            auto n = 10;
            auto vl = layout(comm_, n, Traits::determine());

            TpetraVectord v(vl, 1.);

            TpetraVectord w(vl, 10.);

            TpetraVectord z = max(w, v);

            utopia_test_assert(approxeq(z, w));
        }

        void kokkos_min() {
            auto n = 10;
            auto vl = layout(comm_, n, Traits::determine());

            TpetraVectord v(vl, 1.);

            TpetraVectord w(vl, 10.);

            TpetraVectord z = min(w, v);

            utopia_test_assert(approxeq(z, v));
        }

        void kokkos_sum_reduction() {
            auto n = 10;
            auto vl = layout(comm_, n, Traits::determine());

            TpetraVectord v(vl, 1.);

            double z = sum(v);

            utopia_test_assert(approxeq(z, size(v).get(0) * 1.));
        }

        void kokkos_min_reduction() {
            auto n = 10;
            auto vl = layout(comm_, n, Traits::determine());

            TpetraVectord v(vl, 1.);

            auto r = range(v);

            auto v_view = view_device(v);
            parallel_for(
                range_device(v), UTOPIA_LAMBDA(const SizeType &i) { v_view.set(i, i - r.begin()); });

            Scalar z = min(v);

            utopia_test_assert(approxeq(0., z));
        }

        void kokkos_max_reduction() {
            auto n = 10;
            auto vl = layout(comm_, n, Traits::determine());

            TpetraVectord v(vl, 1.);

            auto r = range(v);

            auto v_view = view_device(v);
            parallel_for(
                range_device(v), UTOPIA_LAMBDA(const SizeType &i) { v_view.set(i, i - r.begin()); });

            Scalar z = max(v);

            utopia_test_assert(approxeq(9., z));
        }

        void kokkos_write() {
            auto n = 10;
            auto vl = layout(comm_, n, Traits::determine());

            TpetraVectord w(vl, -1.);

            auto w_view = view_device(w);
            parallel_for(
                range_device(w), UTOPIA_LAMBDA(const SizeType &i) { w_view.set(i, i); });

            {
                Read<TpetraVectord> r_(w);
                auto r = range(w);

                for (auto i = r.begin(); i < r.end(); ++i) {
                    utopia_test_assert(approxeq(w.get(i), Scalar(i)));
                }
            }
        }

        void kokkos_parallel_each_mat() {
            auto n = 10;
            auto ml = layout(comm_, n, n, Traits::determine(), Traits::determine());

            TpetraMatrixd w;
            w.identity(ml, 1.0);

            w.transform_ijv(
                UTOPIA_LAMBDA(const SizeType &i, const SizeType &j, const Scalar &)->Scalar { return i * n + j; });

            // serial implementation for test
            w.read(UTOPIA_LAMBDA(const SizeType &i, const SizeType &j, const Scalar &val) {
                // We use assert for the gpu
                UTOPIA_DEVICE_ASSERT(device::approxeq(Scalar(i * n + j), val));
            });
        }

        void kokkos_read() {
            auto n = 10;
            auto vl = layout(comm_, n, Traits::determine());

            TpetraVectord w(vl, 50);

            auto w_view = const_view_device(w);
            parallel_for(
                range_device(w), UTOPIA_LAMBDA(const SizeType &i) {
                    auto v = w_view.get(i);
                    UTOPIA_UNUSED(v);
                });
        }

        void kokkos_apply() {
            auto nr = 3;
            auto nc = 3;

            auto ml = layout(comm_, nr, nc, Traits::determine(), Traits::determine());

            TpetraMatrixd P;
            P.sparse(ml, 1, 0);

            {
                Write<TpetraMatrixd> w_(P);
                auto r = row_range(P);
                auto cols = size(P).get(1);
                for (auto i = r.begin(); i < r.end(); ++i) {
                    if (i >= cols) {
                        break;
                    }

                    P.set(i, i, 1.);
                }
            }

            P.transform_values(UTOPIA_LAMBDA(const Scalar &value)->Scalar { return value * 2.; });
        }

        void run() {
            UTOPIA_RUN_TEST(kokkos_max);
            UTOPIA_RUN_TEST(kokkos_min);
            UTOPIA_RUN_TEST(kokkos_min_reduction);
            UTOPIA_RUN_TEST(kokkos_max_reduction);
            UTOPIA_RUN_TEST(kokkos_sum_reduction);
            UTOPIA_RUN_TEST(kokkos_write);
            UTOPIA_RUN_TEST(kokkos_read);
            UTOPIA_RUN_TEST(kokkos_apply);
            UTOPIA_RUN_TEST(kokkos_parallel_each_mat);
        }

    private:
        Comm comm_;
    };

    static void kokkos() { KokkosTest().run(); }

    UTOPIA_REGISTER_TEST_FUNCTION(kokkos);
}  // namespace utopia

#endif  // WITH_TRILINOS
