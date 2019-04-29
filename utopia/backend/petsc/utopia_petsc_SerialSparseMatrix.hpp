
#ifndef UTOPIA_UTOPIA_SERIAL_SPARSE_PETSC_SPARSEMATRIX_HPP
#define UTOPIA_UTOPIA_SERIAL_SPARSE_PETSC_SPARSEMATRIX_HPP

#include "utopia_petsc_SparseMatrix.hpp"
#include "petscmat.h"

namespace utopia {
    class PetscSerialSparseMatrix : public PetscSparseMatrix {
    public:
        virtual MatType type_override() const override
        {
            return MATSEQAIJ;
        }
    };

}

#endif //UTOPIA_UTOPIA_SERIAL_SPARSE_PETSC_SPARSEMATRIX_HPP
