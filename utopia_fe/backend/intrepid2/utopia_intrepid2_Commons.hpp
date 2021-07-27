#ifndef UTOPIA_INTREPID2_COMMONS_HPP
#define UTOPIA_INTREPID2_COMMONS_HPP

#include "utopia_Base.hpp"
#include "utopia_Instance.hpp"

#include <Kokkos_Core.hpp>
#include <Kokkos_View.hpp>

namespace utopia {
    namespace intrepid2 {
        template <class View>
        void fill(View &v, const typename View::const_data_type &value) {
            using ExecutionSpace_t = typename View::execution_space;

            switch (v.rank()) {
                case 1: {
                    Kokkos::parallel_for(
                        v.extent(0), UTOPIA_LAMBDA(const int i) { v(i) = value; });
                    return;
                }
                case 2: {
                    Kokkos::parallel_for(
                        Kokkos::MDRangePolicy<Kokkos::Rank<2>, ExecutionSpace_t>({0, 0}, {v.extent(0), v.extent(1)}),
                        UTOPIA_LAMBDA(const int i, const int j) { v(i, j) = value; });
                    return;
                }

                case 3: {
                    Kokkos::parallel_for(
                        Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace_t>(
                            {0, 0, 0}, {v.extent(0), v.extent(1), v.extent(2)}),
                        UTOPIA_LAMBDA(const int i, const int j, const int k) { v(i, j, k) = value; });
                    return;
                }

                case 4: {
                    Kokkos::parallel_for(
                        Kokkos::MDRangePolicy<Kokkos::Rank<4>, ExecutionSpace_t>(
                            {0, 0, 0, 0}, {v.extent(0), v.extent(1), v.extent(2), v.extent(3)}),
                        UTOPIA_LAMBDA(const int i, const int j, const int k, const int l) { v(i, j, k, l) = value; });
                    return;
                }

                case 5: {
                    Kokkos::parallel_for(
                        Kokkos::MDRangePolicy<Kokkos::Rank<5>, ExecutionSpace_t>(
                            {0, 0, 0, 0, 0}, {v.extent(0), v.extent(1), v.extent(2), v.extent(3), v.extent(4)}),
                        UTOPIA_LAMBDA(const int i, const int j, const int k, const int l, const int m) {
                            v(i, j, k, l, m) = value;
                        });
                    return;
                }
                default: {
                    assert(false);
                    Utopia::Abort("fill: Only up to Rank 5 Views are supported!");
                }
            }
        }

    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_COMMONS_HPP
