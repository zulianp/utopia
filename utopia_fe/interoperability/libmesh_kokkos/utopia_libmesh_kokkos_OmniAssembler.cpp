#include "utopia_libmesh_kokkos_OmniAssembler.hpp"
#include "utopia_kokkos_OmniAssembler_impl.hpp"
#include "utopia_libmesh_kokkos.hpp"

namespace utopia {
    namespace kokkos {
        template class OmniAssembler<utopia::libmesh::FunctionSpace, utopia::kokkos::FE<double>>;
    }
}  // namespace utopia
