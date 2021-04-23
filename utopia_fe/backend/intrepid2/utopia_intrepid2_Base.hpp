#ifndef UTOPIA_INTREPID2_BASE_HPP
#define UTOPIA_INTREPID2_BASE_HPP

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

        using IntViewDevice = utopia::intrepid2::ViewDevice<int>;
        using IntViewHost = utopia::intrepid2::ViewHost<int>;
    }  // namespace intrepid2

    template <typename Scalar_, typename... Args>
    class Traits<intrepid2::View<Scalar_, Args...>> {
    public:
        using Scalar = Scalar_;
    };

    // template <typename Scalar_>
    // class Traits<intrepid2::ViewHost<Scalar_>> {
    // public:
    //     using Scalar = Scalar_;
    // };

}  // namespace utopia

#endif