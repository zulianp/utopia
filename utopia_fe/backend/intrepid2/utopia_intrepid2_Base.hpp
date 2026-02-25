#ifndef UTOPIA_INTREPID2_BASE_HPP
#define UTOPIA_INTREPID2_BASE_HPP

#include "utopia_Traits.hpp"

#include <Kokkos_DynRankView.hpp>

namespace utopia {
    namespace intrepid2 {
        using HostExecutionSpace = ::Kokkos::DefaultHostExecutionSpace;
        using ExecutionSpace = ::Kokkos::DefaultExecutionSpace;

        template <typename T, typename... Args>
        using View = ::Kokkos::DynRankView<T, Args...>;

        template <typename T>
        using ViewDevice = ::Kokkos::DynRankView<T, ExecutionSpace>;

        template <typename T>
        using ViewHost = ::Kokkos::DynRankView<T, HostExecutionSpace>;

#if KOKKOS_VERSION >= 50002
        template <typename T>
        using ViewHostMirror = typename ViewDevice<T>::host_mirror_type;
        using IntViewHostMirror = utopia::intrepid2::ViewDevice<int>::host_mirror_type;

#else
        template <typename T>
        using ViewHostMirror = typename ViewDevice<T>::HostMirror;
        using IntViewHostMirror = utopia::intrepid2::ViewDevice<int>::HostMirror;
#endif


        using IntViewDevice = utopia::intrepid2::ViewDevice<int>;
        using IntViewHost = utopia::intrepid2::ViewHost<int>;
    }  // namespace intrepid2

    // template <typename Scalar_, typename... Args>
    // class Traits<intrepid2::View<Scalar_, Args...>> {
    // public:
    //     using Scalar = Scalar_;
    // };

    // template <typename Scalar_>
    // class Traits<intrepid2::ViewHost<Scalar_>> {
    // public:
    //     using Scalar = Scalar_;
    // };

}  // namespace utopia

#endif