#include "utopia_petsc_SideSet.hpp"

namespace utopia {
    namespace petsc {
        constexpr SideSets<1>::Sides SideSets<1>::ids_;
        constexpr SideSets<2>::Sides SideSets<2>::ids_;
        constexpr SideSets<3>::Sides SideSets<3>::ids_;
    }  // namespace petsc
}  // namespace utopia
