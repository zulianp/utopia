#ifndef UTOPIA_KOKKOS_TRAITS_HPP
#define UTOPIA_KOKKOS_TRAITS_HPP

#include "utopia_Algorithms.hpp"
#include "utopia_BackendInfo.hpp"
#include "utopia_MatrixView.hpp"
#include "utopia_Traits.hpp"
#include "utopia_VectorView.hpp"
#include "utopia_kokkos_ForwardDeclarations.hpp"
#include "utopia_trilinos_Base.hpp"

#include <KokkosBlas1_axpby.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosBlas1_fill.hpp>
#include <KokkosBlas1_scal.hpp>
#include <KokkosBlas2_gemv.hpp>
#include <KokkosBlas3_gemm.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_View.hpp>

namespace utopia {

    namespace device {

        template <typename Scalar, class... Args>
        class Fill<Kokkos::View<Scalar **, Args...>> {
        public:
            using View = Kokkos::View<Scalar **, Args...>;

            UTOPIA_INLINE_FUNCTION static void apply(const Scalar &val, View &in_out) {
                const SizeType rows = in_out.extent(0);
                const SizeType cols = in_out.extent(1);

                for (SizeType i = 0; i < rows; ++i) {
                    for (SizeType j = 0; j < cols; ++j) {
                        in_out(i, j) = val;
                    }
                }
            }
        };

        template <class Scalar, class... Args>
        class AXPY<Kokkos::View<Scalar **, Args...>, Kokkos::View<Scalar **, Args...>> {
        public:
            using View = Kokkos::View<Scalar **, Args...>;

            UTOPIA_INLINE_FUNCTION static void apply(const Scalar &alpha, const View &x, View &y) {
                KokkosBlas::axpy(alpha, x, y);
            }
        };

        template <class Scalar, class... Args>
        class Copy<Kokkos::View<Scalar **, Args...>, Kokkos::View<Scalar **, Args...>> {
        public:
            using View = Kokkos::View<Scalar **, Args...>;

            UTOPIA_INLINE_FUNCTION static void apply(const View &in, View &out) {
                const SizeType rows = in.extent(0);
                const SizeType cols = in.extent(1);

                for (SizeType i = 0; i < rows; ++i) {
                    for (SizeType j = 0; j < cols; ++j) {
                        out(i, j) = in(i, j);
                    }
                }
            }
        };

    }  // namespace device

    template <class Scalar, typename... Args>
    class Accessor<Kokkos::View<Scalar **, Args...>> {
    public:
        using View = Kokkos::View<Scalar **, Args...>;
        using SizeType = utopia::TpetraSizeType;

        UTOPIA_INLINE_FUNCTION static const Scalar &get(const View &t, const SizeType &i, const SizeType &j) {
            return t(i, j);
        }

        UTOPIA_INLINE_FUNCTION static void set(View &t, const SizeType &i, const SizeType &j, const Scalar &val) {
            t(i, j) = val;
        }
    };

    template <class Scalar, typename... Args>
    class Accessor<Kokkos::View<Scalar *, Args...>> {
    public:
        using View = Kokkos::View<Scalar *, Args...>;
        using SizeType = utopia::TpetraSizeType;

        UTOPIA_INLINE_FUNCTION static const Scalar &get(const View &t, const SizeType &i) { return t(i); }

        UTOPIA_INLINE_FUNCTION static void set(View &t, const SizeType &i, const Scalar &val) { t(i) = val; }
    };

    // template<
    //     class Matrix_,
    //     class Vector_
    //     class Scalar_,
    //     class SizeType_>
    // class KokkosBaseTraits {
    // public:
    //     using Scalar   = Scalar_;
    //     using SizeType = SizeType_;
    //     using Vector   = Vector_;
    //     using Matrix   = Matrix_;

    //     //FIXME use Kokkos compatible wrapper
    //     using IndexSet    = utopia::TpetraIndexSet;
    //     using IndexArray  = utopia::TpetraIndexArray;
    //     using ScalarArray = utopia::TpetraScalarArray;

    //     enum {
    //         Backend = KOKKOS
    //     };

    //     static BackendInfo &backend_info()
    //     {
    //         static BackendInfo instance_("kokkos");
    //         return instance_;
    //     }
    // };

    template <class Scalar_>
    class KokkosTraits {
    public:
        using Scalar = Scalar_;
        using SizeType = utopia::TpetraSizeType;

        using KokkosView1 = Kokkos::View<Scalar *>;
        using Vector = utopia::VectorView<KokkosView1>;

        using KokkosView2 = Kokkos::View<Scalar **>;
        using Matrix = utopia::MatrixView<KokkosView2>;

        // FIXME use Kokkos compatible wrapper
        using IndexSet = utopia::TpetraIndexSet;
        using IndexArray = utopia::TpetraIndexArray;
        using ScalarArray = utopia::TpetraScalarArray;

        enum { Backend = KOKKOS };

        static BackendInfo &backend_info() {
            static BackendInfo instance_("kokkos");
            return instance_;
        }
    };

    template <class Scalar_, typename... Args>
    class Traits<Kokkos::View<Scalar_ *, Args...>> {
    public:
        using Scalar = Scalar_;
        using SizeType = std::size_t;
    };

    template <class Scalar_, typename... Args>
    class Traits<Kokkos::View<Scalar_ **, Args...>> {
    public:
        using Scalar = Scalar_;
        using SizeType = std::size_t;
    };

    template <typename Scalar>
    using DefaultVectorView = utopia::VectorView<Kokkos::View<Scalar *>>;

    template <typename Scalar>
    using DefaultMatrixView = utopia::MatrixView<Kokkos::View<Scalar **>>;

    UTOPIA_MAKE_TRAITS_TPL_1(DefaultVectorView, KokkosTraits, 1);
    UTOPIA_MAKE_TRAITS_DENSE_TPL_1(DefaultMatrixView, KokkosTraits, 2);
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_TRAITS_HPP
