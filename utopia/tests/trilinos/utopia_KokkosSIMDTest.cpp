#include "utopia_Base.hpp"

#ifdef UTOPIA_WITH_KOKKOS_SIMD

#include "utopia.hpp"
#include "utopia_Testing.hpp"

#include "simd.hpp"
#include "utopia_Algorithms.hpp"
#include "utopia_Utils.hpp"

#include "utopia_Algorithms.hpp"
#include "utopia_Trace.hpp"
#include "utopia_VectorView.hpp"
#include "utopia_Views.hpp"

#include <type_traits>
#include <vector>

namespace utopia {

    template <class Vector, int Backend = Traits<Vector>::Backend>
    class JacobiSweep {};

    template <class Vector>
    class JacobiSweep<Vector, TRILINOS> {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using LocalSizeType = typename Traits<Vector>::LocalSizeType;
        using ViewType = Kokkos::View<Scalar *>;

        using SIMDType = simd::native_simd<Scalar>;
        static const int SIMD_LANES = SIMDType::size();

        using SIMDViewType = Kokkos::View<SIMDType *>;

        void init(const Vector &t_inv_diag) {
            const auto &inv_diag = t_inv_diag.derived();
            const SizeType n_rows = inv_diag.local_size();
            const SizeType n_blocks = n_rows / SIMD_LANES;
            const SizeType n_clean = n_blocks * SIMD_LANES;
            const SizeType rmd = n_rows - n_clean;
            const SizeType n_blocks_with_rmd = n_blocks + (rmd > 0);

            if (n_blocks_with_rmd != SizeType(d_inv_.size())) {
                d_inv_ = SIMDViewType("d_inv", n_blocks_with_rmd);
                c1_ = SIMDViewType("c1", n_blocks_with_rmd);
                c2_ = SIMDViewType("c2", n_blocks_with_rmd);
            }

            copy(inv_diag, d_inv_);
        }

        void apply(const Vector &in, Vector &out) {
            // Vectorize input
            copy(in, c1_);

            const SizeType n_rows = d_inv_.size();

            for (SizeType t = 0; t < n_sweeps_; ++t) {
                if (t % 2 == 0) {
                    Kokkos::parallel_for(
                        n_rows, KOKKOS_LAMBDA(const int i) { c2_[i] = d_inv_[i] * c1_[i]; });
                } else {
                    Kokkos::parallel_for(
                        n_rows, KOKKOS_LAMBDA(const int i) { c1_[i] = d_inv_[i] * c2_[i]; });
                }
            }

            // Copy back
            if (n_sweeps_ % 2 == 0) {
                copy(c1_, out);
            } else {
                copy(c2_, out);
            }
        }

        JacobiSweep() = default;

        inline void sweeps(const SizeType n) { n_sweeps_ = n; }

    private:
        SIMDViewType c1_;
        SIMDViewType c2_;
        SIMDViewType d_inv_;
        // FIXME
        SizeType n_sweeps_{100};

        static void copy(const Vector &in, SIMDViewType &out) {
            const SizeType n = in.local_size();
            const SizeType n_blocks = n / SIMD_LANES;
            const SizeType n_clean = n_blocks * SIMD_LANES;
            const SizeType rmd = n - n_clean;

            auto view = local_view_device(in);

            Kokkos::parallel_for(
                n_blocks, KOKKOS_LAMBDA(const int i) {
                    Scalar buff[SIMD_LANES];

                    LocalSizeType idx = i * SIMD_LANES;

                    for (int k = 0; k < SIMD_LANES; ++k) {
                        buff[k] = view.get(idx + k);
                    }

                    out(i).copy_from(buff, simd::element_aligned_tag());
                });

            if (rmd) {
                Kokkos::parallel_for(
                    1, KOKKOS_LAMBDA(const int &) {
                        Scalar buff[SIMD_LANES];

                        LocalSizeType idx = n_blocks * SIMD_LANES;

                        for (int k = 0; k < rmd; ++k) {
                            buff[k] = view.get(idx + k);
                        }

                        for (int k = rmd; k < SIMD_LANES; ++k) {
                            buff[k] = 0;
                        }

                        out(n_blocks).copy_from(buff, simd::element_aligned_tag());
                    });
            }
        }

        static void copy(const SIMDViewType &in, Vector &out) {
            const SizeType n = out.local_size();
            const SizeType n_blocks = n / SIMD_LANES;
            const SizeType n_clean = n_blocks * SIMD_LANES;
            const SizeType rmd = n - n_clean;

            auto view = local_view_device(out);

            Kokkos::parallel_for(
                n_blocks, KOKKOS_LAMBDA(const int i) {
                    Scalar buff[SIMD_LANES];

                    LocalSizeType idx = i * SIMD_LANES;

                    in(i).copy_to(buff, simd::element_aligned_tag());

                    for (int k = 0; k < SIMD_LANES; ++k) {
                        view.set(idx + k, buff[k]);
                    }
                });

            if (rmd) {
                Kokkos::parallel_for(
                    1, KOKKOS_LAMBDA(const int &) {
                        Scalar buff[SIMD_LANES];

                        LocalSizeType idx = n_blocks * SIMD_LANES;

                        in(n_blocks).copy_to(buff, simd::element_aligned_tag());

                        for (int k = 0; k < rmd; ++k) {
                            view.set(idx + k, buff[k]);
                        }
                    });
            }
        }
    };

    template class JacobiSweep<TpetraVector, TRILINOS>;
}  // namespace utopia

using namespace utopia;

static void kokkos_simd_jacobi() {
    using Scalar = Traits<TpetraVector>::Scalar;
    using SizeType = Traits<TpetraVector>::SizeType;

    auto &comm = TrilinosCommunicator::get_default();
    SizeType n = 10000000;
    auto vl = layout(comm, n, n * comm.size());

    JacobiSweep<TpetraVector> jac;
    TpetraVector in, out, inv_diag;

    auto sweeps = 100;
    jac.sweeps(sweeps);

    in.values(vl, std::pow(2.0, sweeps));
    out.values(vl, 0.0);
    inv_diag.values(vl, 0.5);

    Chrono c;
    c.start();
    jac.init(inv_diag);
    c.stop();

    if (comm.rank() == 0) std::cout << c << std::endl;

    c.start();
    jac.apply(in, out);
    c.stop();

    if (comm.rank() == 0) std::cout << c << std::endl;
    // disp(out);

    Scalar actual = sum(out);

    disp(actual);

    utopia_test_assert(approxeq(actual, static_cast<Scalar>(out.size()), device::epsilon<Scalar>()));
}

UTOPIA_REGISTER_TEST_FUNCTION(kokkos_simd_jacobi);

static void kokkos_simd() {
    static const int N = simd::native_simd<double>::size();
    using utopia::device::approxeq;

    double arr_x[N];
    double arr_y[N];

    for (int i = 0; i < N; ++i) {
        arr_x[i] = 1;
        arr_y[i] = 2;
    }

    simd::native_simd<double> x(arr_x, simd::element_aligned_tag());
    simd::native_simd<double> y(arr_y, simd::element_aligned_tag());

    y += 0.1 * x;

    y.copy_to(arr_y, simd::element_aligned_tag());

    for (int i = 0; i < N; ++i) {
        // std::cout << arr_y[i] << std::endl;
        utopia_test_assert(approxeq(arr_y[i], 2.1, device::epsilon<double>()));
    }

    disp("simd::native_simd<double>");
    disp(simd::native_simd<double>::size());

    disp("simd::native_simd<float>");
    disp(simd::native_simd<float>::size());
}

UTOPIA_REGISTER_TEST_FUNCTION(kokkos_simd);

#endif  // UTOPIA_WITH_KOKKOS_SIMD