#ifndef UTOPIA_PETSC_INTREPID2_OMNI_ASSEMBLER_HPP
#define UTOPIA_PETSC_INTREPID2_OMNI_ASSEMBLER_HPP

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_OmniAssembler.hpp"
#include "utopia_petsc_FunctionSpace.hpp"

#include "utopia_kokkos_OmniAssembler_impl.hpp"

#include "utopia_petsc_intrepid2.hpp"

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_OmniAssembler.hpp"
#include "utopia_petsc_FunctionSpace.hpp"

namespace utopia {

    template <>
    class OmniAssembler<utopia::petsc::FunctionSpace> : public intrepid2::OmniAssembler<utopia::petsc::FunctionSpace> {
    public:
        using Super = intrepid2::OmniAssembler<utopia::petsc::FunctionSpace>;
        using Super::Super;
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_INTREPID2_OMNI_ASSEMBLER_HPP