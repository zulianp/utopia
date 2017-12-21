
#ifndef UTOPIA_PETSC_SPARSE_MATRIX_HPP
#define UTOPIA_PETSC_SPARSE_MATRIX_HPP

#include "utopia_petsc_Matrix.hpp"
#include <map>

namespace utopia{

	class PETScSparseMatrix : public PETScMatrix {
	public:
		using PETScMatrix::PETScMatrix;
		using PETScMatrix::operator=;
		virtual ~PETScSparseMatrix() {}
	};

}

#endif //UTOPIA_PETSC_SPARSE_MATRIX_HPP
