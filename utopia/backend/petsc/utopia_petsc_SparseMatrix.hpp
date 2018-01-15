
#ifndef UTOPIA_PETSC_SPARSE_MATRIX_HPP
#define UTOPIA_PETSC_SPARSE_MATRIX_HPP

#include "utopia_petsc_Matrix.hpp"
#include <map>

namespace utopia{

	class PetscSparseMatrix : public PetscMatrix {
	public:
		using PetscMatrix::PetscMatrix;
		using PetscMatrix::operator=;
		virtual ~PetscSparseMatrix() {}
	};

}

#endif //UTOPIA_PETSC_SPARSE_MATRIX_HPP
