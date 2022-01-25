#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_UniformFE.hpp"

namespace utopia {
    namespace kokkos {
        template class UniformFE<double>;
        template class FE<double>;
    }  // namespace kokkos
}  // namespace utopia
