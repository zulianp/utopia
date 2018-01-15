
#ifndef UTOPIA_UTOPIA_SERIAL_SPARSE_PETSC_SPARSEMATRIX_HPP
#define UTOPIA_UTOPIA_SERIAL_SPARSE_PETSC_SPARSEMATRIX_HPP

#include "utopia_petsc_SparseMatrix.hpp"

namespace utopia {
	class PetscSerialSparseMatrix : public PetscSparseMatrix {
	public:
		using PetscSparseMatrix::PetscSparseMatrix;
		using PetscSparseMatrix::operator=;
	};

}

#endif //UTOPIA_UTOPIA_SERIAL_SPARSE_PETSC_SPARSEMATRIX_HPP
