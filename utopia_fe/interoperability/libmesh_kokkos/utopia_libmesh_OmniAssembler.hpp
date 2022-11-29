#ifndef UTOPIA_PETSC_INTREPID2_OMNI_ASSEMBLER_HPP
#define UTOPIA_PETSC_INTREPID2_OMNI_ASSEMBLER_HPP

#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_OmniAssembler.hpp"
#include "utopia_libmesh_FunctionSpace_new.hpp"

#include "utopia_kokkos_OmniAssembler_impl.hpp"

#include "utopia_libmesh_kokkos.hpp"

#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_OmniAssembler.hpp"

namespace utopia {

    template <>
    class OmniAssembler<utopia::libmesh::FunctionSpace>
        : public kokkos::OmniAssembler<utopia::libmesh::FunctionSpace,
                                       utopia::kokkos::FE<Traits<utopia::libmesh::FunctionSpace>::Scalar>> {
    public:
        using Scalar = Traits<utopia::libmesh::FunctionSpace>::Scalar;
        using FE = utopia::kokkos::FE<Scalar>;
        using Super = kokkos::OmniAssembler<utopia::libmesh::FunctionSpace, FE>;
        using Super::Super;
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_INTREPID2_OMNI_ASSEMBLER_HPP