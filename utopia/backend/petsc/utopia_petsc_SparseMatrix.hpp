
#ifndef UTOPIA_PETSC_SPARSE_MATRIX_HPP
#define UTOPIA_PETSC_SPARSE_MATRIX_HPP

#include "utopia_petsc_Matrix.hpp"
#include <map>

namespace utopia{

	class PETScSparseMatrix : public PETScMatrix {
	public:
		virtual ~PETScSparseMatrix() {}

		//FIXME is this really necessary???
		std::map<std::pair<PetscInt, PetscInt>, PetscScalar> buffer;
	};

}

#endif //UTOPIA_PETSC_SPARSE_MATRIX_HPP
