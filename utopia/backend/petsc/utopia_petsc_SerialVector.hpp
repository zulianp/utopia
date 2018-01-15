#ifndef UTOPIA_UTOPIA_SERIAL_SPARSE_PETSC_SPARSEVECTOR_HPP
#define UTOPIA_UTOPIA_SERIAL_SPARSE_PETSC_SPARSEVECTOR_HPP

#include "utopia_petsc_Vector.hpp"

namespace utopia {
	class PetscSerialVector : public PetscVector {
	public:
		using PetscVector::PetscVector;
		using PetscVector::operator=;
	};
}

#endif //UTOPIA_UTOPIA_SERIAL_SPARSE_PETSC_SPARSEVECTOR_HPP
