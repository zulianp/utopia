#ifndef UTOPIA_KOKKOS_PARALLEL_FOR_HPP
#define UTOPIA_KOKKOS_PARALLEL_FOR_HPP

#include "utopia_Base.hpp"

#include "utopia_For.hpp"
#include "utopia_Traits.hpp"

#include <Kokkos_Core.hpp>

namespace utopia {

    template <>
    class ParallelFor<TRILINOS> {
    public:
        template <typename F>
        inline static void apply(const Range &r, F f) {
            apply(r.begin(), r.end(), f);
        }

        template <typename F>
        inline static void apply(const std::size_t &begin, const std::size_t &end, F f) {
            auto extent = end - begin;
            Kokkos::parallel_for(
                extent, KOKKOS_LAMBDA(const int i) { f(begin + i); });
        }

        template <typename F>
        inline static void apply(const std::size_t &n, F f) {
            Kokkos::parallel_for(
                n, KOKKOS_LAMBDA(const int i) { f(i); });
        }
    };
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_PARALLEL_FOR_HPP
