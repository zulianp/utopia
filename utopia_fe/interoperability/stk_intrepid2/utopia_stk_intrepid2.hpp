#ifndef UTOPIA_STK_INTREPID2_HPP
#define UTOPIA_STK_INTREPID2_HPP

#include "utopia_Base.hpp"

#ifdef UTOPIA_WITH_PETSC
#include "utopia_petsc_ForwardDeclarations.hpp"
#endif

#include "utopia_CreateFE.hpp"
#include "utopia_LocalToGlobal.hpp"

#include "utopia_intrepid2_ForwardDeclarations.hpp"
#include "utopia_stk_ForwardDeclarations.hpp"

#include <Kokkos_DynRankView.hpp>

namespace utopia {

    template <typename Scalar>
    class CreateFE<utopia::stk::FunctionSpace, utopia::intrepid2::FE<Scalar>> {
    public:
        static void apply(const utopia::stk::FunctionSpace &space,
                          utopia::intrepid2::FE<Scalar> &fe,
                          const int degree = 0);
    };

#ifdef UTOPIA_WITH_PETSC
    template <typename Scalar>
    class LocalToGlobal<utopia::stk::FunctionSpace, ::Kokkos::DynRankView<Scalar>, PetscMatrix> {
    public:
        using DynRankView = ::Kokkos::DynRankView<Scalar>;
        static void apply(const utopia::stk::FunctionSpace &space,
                          const DynRankView &element_matrices,
                          AssemblyMode mode,
                          PetscMatrix &matrix);
    };

    template <typename Scalar>
    class LocalToGlobal<utopia::stk::FunctionSpace, ::Kokkos::DynRankView<Scalar>, PetscVector> {
    public:
        using DynRankView = ::Kokkos::DynRankView<Scalar>;
        static void apply(const utopia::stk::FunctionSpace &space,
                          const DynRankView &element_vectors,
                          AssemblyMode mode,
                          PetscVector &vector);
    };
#endif  // UTOPIA_WITH_PETSC

}  // namespace utopia

#endif  // UTOPIA_STK_INTREPID2_HPP