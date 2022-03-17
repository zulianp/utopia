#include "utopia_assembly_kokkos_generated.hpp"

#include "utopia_fe_kokkos_generated.hpp"
#include "utopia_material_Mass_2_impl.hpp"

namespace utopia {
    namespace kernels {
        // explicit instantiations declaratins
        template class Mass<Quad4<double, double>>;
    }  // namespace kernels
}  // namespace utopia
