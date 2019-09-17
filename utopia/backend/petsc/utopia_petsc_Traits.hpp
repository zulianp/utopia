
#ifndef UTOPIA_UTOPIA_PETSC_TRAITS_HPP
#define UTOPIA_UTOPIA_PETSC_TRAITS_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_BackendInfo.hpp"

#include "utopia_petsc_Base.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"

#include "petscsys.h"
#include "petscmat.h"

namespace utopia {


    class PetscTraits {
    public:
        using Scalar   = PetscScalar;
        using SizeType = PetscInt;

        using Matrix            = utopia::PetscMatrix;
        using SparseMatrix      = utopia::PetscMatrix;
        using PolymorphicMatrix = utopia::PetscMatrix;
        using Vector            = utopia::PetscVector;

        using IndexSet    = utopia::PetscIndexSet;
        using IndexArray  = utopia::PetscArray<SizeType>;
        using ScalarArray = utopia::PetscArray<Scalar>;

        enum
        {
            Backend = PETSC
        };

        static BackendInfo &backend_info()
        {
            static BackendInfo instance_("petsc");
            return instance_;
        }
    };

    UTOPIA_MAKE_TRAITS_POLYMORPHIC(PetscMatrix, PetscTraits);
    UTOPIA_MAKE_TRAITS(PetscVector, PetscTraits);
    
    UTOPIA_MAKE_PARALLEL_TRAITS(PetscMatrix);
    UTOPIA_MAKE_PARALLEL_TRAITS(PetscVector);
}

#endif //UTOPIA_UTOPIA_PETSCTRAITS_HPP
