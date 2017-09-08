
#ifndef UTOPIA_PETSC_SPARSE_MATRIX_HPP
#define UTOPIA_PETSC_SPARSE_MATRIX_HPP

#include "utopia_PETScMatrix.hpp"

namespace utopia{

	class PETScSparseMatrix : public PETScMatrix {
	public:
		virtual ~PETScSparseMatrix() {}

	};

}

#endif //UTOPIA_PETSC_SPARSE_MATRIX_HPP
