#ifndef UTOPIA_KOKKOS_PARALLEL_REDUCE_HPP
#define UTOPIA_KOKKOS_PARALLEL_REDUCE_HPP

#include "utopia_Base.hpp"

#include "utopia_Reduction.hpp"
#include "utopia_Traits.hpp"

#include <Kokkos_Core.hpp>

namespace utopia {

    template <>
    class ParallelReduce<TRILINOS> {
    public:
        template <typename F, class Accumulator>
        inline static void apply(const Range &r, F f, Accumulator &acc) {
            apply(r.begin(), r.end(), f, acc);
        }

        template <typename F, class Accumulator>
        inline static void apply(const std::size_t &begin, const std::size_t &end, F f, Accumulator &acc) {
            apply(end - begin, f, acc);
        }

        template <typename F, class Accumulator>
        inline static void apply(const std::size_t &n, F f, Accumulator &acc) {
            Kokkos::parallel_reduce(
                n, KOKKOS_LAMBDA(const int &i, double &lsum) { lsum += f(i); }, acc);

            Kokkos::fence();
        }
    };
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_PARALLEL_REDUCE_HPP
