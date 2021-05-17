#include "utopia_libmesh_intrepid2_OmniAssembler.hpp"
#include "utopia_intrepid2_OmniAssembler_impl.hpp"
#include "utopia_libmesh_intrepid2.hpp"

namespace utopia {
    namespace intrepid2 {
        template class OmniAssembler<utopia::libmesh::FunctionSpace>;
    }
}  // namespace utopia
