
#ifndef UTOPIA_UTOPIA_PETSCTRAITS_HPP
#define UTOPIA_UTOPIA_PETSCTRAITS_HPP

#include "utopia_Traits.hpp"

#include "utopia_PETScMatrix.hpp"
#include "utopia_PETScSparseMatrix.hpp"
#include "utopia_PETScSerialSparseMatrix.hpp"

#include "utopia_PETScVector.hpp"
#include "utopia_PETScSerialVector.hpp"

#include "utopia_Base.hpp"

#include "petscsys.h"
#include "petscmat.h"

namespace utopia 
{
    class PETScTraits 
    {
    public:
        typedef PetscScalar Scalar;
        typedef PetscInt SizeType;

        typedef utopia::PETScMatrix Matrix;
        typedef utopia::PETScSparseMatrix SparseMatrix;
        typedef utopia::PETScSerialSparseMatrix SerialSparseMatrix;

        typedef utopia::PETScVector Vector;
        typedef utopia::PETScSerialVector SerialVector;

        enum 
        {
            Backend = PETSC
        };
    };
}

#endif //UTOPIA_UTOPIA_PETSCTRAITS_HPP
