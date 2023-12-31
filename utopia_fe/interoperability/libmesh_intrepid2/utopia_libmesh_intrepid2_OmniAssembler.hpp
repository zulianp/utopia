#ifndef UTOPIA_LIBMESH_INTREPID2_OMNI_ASSEMBLER_HPP
#define UTOPIA_LIBMESH_INTREPID2_OMNI_ASSEMBLER_HPP

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_fe_Core.hpp"

#include "utopia_fe_Environment.hpp"

#include "utopia_intrepid2_OmniAssembler.hpp"
#include "utopia_libmesh_ForwardDeclarations.hpp"
#include "utopia_libmesh_FunctionSpace_new.hpp"

namespace utopia {

    template <>
    class OmniAssembler<utopia::libmesh::FunctionSpace> final
        : public utopia::intrepid2::OmniAssembler<utopia::libmesh::FunctionSpace> {
    public:
        using Super = utopia::intrepid2::OmniAssembler<utopia::libmesh::FunctionSpace>;
        using Super::Super;
    };

}  // namespace utopia

#endif  // UTOPIA_LIBMESH_INTREPID2_OMNI_ASSEMBLER_HPP
