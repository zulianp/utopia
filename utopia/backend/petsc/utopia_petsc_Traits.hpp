
#ifndef UTOPIA_UTOPIA_PETSCTRAITS_HPP
#define UTOPIA_UTOPIA_PETSCTRAITS_HPP

#include "utopia_Traits.hpp"

#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_SparseMatrix.hpp"
#include "utopia_petsc_SerialSparseMatrix.hpp"

#include "utopia_petsc_Vector.hpp"
#include "utopia_petsc_SerialVector.hpp"

#include "utopia_Base.hpp"

#include "petscsys.h"
#include "petscmat.h"

namespace utopia {
    class PetscTraits {
    public:
        typedef PetscScalar Scalar;
        typedef PetscInt SizeType;

        typedef utopia::PetscMatrix Matrix;
        typedef utopia::PetscSparseMatrix SparseMatrix;
        typedef utopia::PetscSerialSparseMatrix SerialSparseMatrix;

        typedef utopia::PetscVector Vector;
        typedef utopia::PetscSerialVector SerialVector;

        enum
        {
            Backend = PETSC
        };
    };


    class PetscCudaTraits {
    public:
        typedef PetscScalar Scalar;
        typedef PetscInt SizeType;

        typedef utopia::PetscCuSparseMatrix Matrix;
        typedef utopia::PetscCuSparseMatrix SparseMatrix;

        typedef utopia::PetscCuVector Vector;
        typedef utopia::PetscCuVector SerialVector;

        enum
        {
            Backend = PETSC
        };
    };
}

#endif //UTOPIA_UTOPIA_PETSCTRAITS_HPP
