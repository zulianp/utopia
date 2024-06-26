

#include "utopia_fe_kokkos_generated.hpp"

namespace utopia {
    namespace kernels {
        // explicit instantiations declaratins

        template class AxisAlignedHex8<double, double>;

        template class AxisAlignedQuad4<double, double>;

        template class Hex8<double, double>;

        template class Pentatope5<double, double>;

        template class Quad4<double, double>;

        template class Tet4<double, double>;

        template class Tri3<double, double>;
    }  // namespace kernels
}  // namespace utopia
