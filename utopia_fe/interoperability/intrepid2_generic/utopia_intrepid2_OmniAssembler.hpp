#ifndef UTOPIA_INTREPID2_OMNI_ASSEMBLER_HPP
#define UTOPIA_INTREPID2_OMNI_ASSEMBLER_HPP

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_FEAssembler.hpp"
#include "utopia_kokkos_FEAssembler.hpp"
#include "utopia_kokkos_OmniAssembler.hpp"

#include "utopia_fe_Core.hpp"
#include "utopia_fe_Environment.hpp"
#include "utopia_intrepid2_ForwardDeclarations.hpp"

namespace utopia {
    namespace intrepid2 {

        template <class FunctionSpace>
        using OmniAssembler =
            utopia::kokkos::OmniAssembler<FunctionSpace, utopia::intrepid2::FE<typename Traits<FunctionSpace>::Scalar>>;

    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_OMNI_ASSEMBLER_HPP
