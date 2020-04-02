
#ifndef UTOPIA_UTOPIA_PETSC_TRAITS_HPP
#define UTOPIA_UTOPIA_PETSC_TRAITS_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_BackendInfo.hpp"

#include "utopia_petsc_Base.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_Device.hpp"
#include "utopia_Layout.hpp"

#include "petscsys.h"
#include "petscmat.h"

namespace utopia {


    class PetscTraits {
    public:
        using Scalar        = PetscScalar;
        using SizeType      = PetscInt;
        using LocalSizeType = SizeType;

        using Matrix            = utopia::PetscMatrix;
        using SparseMatrix      = utopia::PetscMatrix;
        using PolymorphicMatrix = utopia::PetscMatrix;
        using Vector            = utopia::PetscVector;

        using IndexSet    = utopia::PetscIndexSet;
        using IndexArray  = utopia::PetscArray<SizeType>;
        using ScalarArray = utopia::PetscArray<Scalar>;

        using Communicator = utopia::PetscCommunicator;
        using Device       = utopia::Device<PETSC>;
        using Layout       = utopia::Layout<PetscCommunicator, 1, SizeType>;
        using MatrixLayout = utopia::Layout<PetscCommunicator, 2, SizeType>;

        enum
        {
            Backend = PETSC
        };

        static BackendInfo &backend_info()
        {
            static BackendInfo instance_("petsc");
            return instance_;
        }

        //tells petsc to device local size automatically
        static constexpr SizeType decide()
        {
            return PETSC_DECIDE;
        }

        //tells petsc to compute global size automatically
        static constexpr SizeType determine()
        {
            return PETSC_DETERMINE;
        }
    };

    UTOPIA_MAKE_TRAITS_POLYMORPHIC(PetscMatrix, PetscTraits, 2);
    UTOPIA_MAKE_TRAITS(PetscVector, PetscTraits, 1);

    UTOPIA_MAKE_PARALLEL_TRAITS(PetscMatrix);
    UTOPIA_MAKE_PARALLEL_TRAITS(PetscVector);
}

#endif //UTOPIA_UTOPIA_PETSCTRAITS_HPP
