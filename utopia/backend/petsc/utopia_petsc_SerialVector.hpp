#ifndef UTOPIA_UTOPIA_SERIAL_SPARSE_PETSC_SPARSEVECTOR_HPP
#define UTOPIA_UTOPIA_SERIAL_SPARSE_PETSC_SPARSEVECTOR_HPP

#include "utopia_petsc_Vector.hpp"

namespace utopia {
    class PetscSerialVector : public PetscVector {
    public:
        virtual VecType type_override() const
        {
            return VECSEQ;
        }
    };

    class PetscCuVector : public PetscVector {
    public:
        virtual VecType type_override() const
        {
            return VECCUDA;
        }
    };
}

#endif //UTOPIA_UTOPIA_SERIAL_SPARSE_PETSC_SPARSEVECTOR_HPP
