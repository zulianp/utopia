#ifndef UTOPIA_STK_INTREPID2_OMNI_ASSEMBLER_HPP
#define UTOPIA_STK_INTREPID2_OMNI_ASSEMBLER_HPP

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_fe_Core.hpp"

#include "utopia_fe_Environment.hpp"

#include "utopia_intrepid2_OmniAssembler.hpp"
#include "utopia_stk_ForwardDeclarations.hpp"
#include "utopia_stk_FunctionSpace.hpp"

namespace utopia {

    template <>
    class OmniAssembler<utopia::stk::FunctionSpace> final
        : public utopia::intrepid2::OmniAssembler<utopia::stk::FunctionSpace> {
    public:
        using Super = utopia::intrepid2::OmniAssembler<utopia::stk::FunctionSpace>;
        using Super::Super;
    };

}  // namespace utopia

#endif  // UTOPIA_STK_INTREPID2_OMNI_ASSEMBLER_HPP
