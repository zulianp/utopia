#ifndef UTOPIA_PETSC_FACTORIZATIONS_HPP
#define UTOPIA_PETSC_FACTORIZATIONS_HPP 

#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_petsc_Factorization.hpp"

namespace utopia {

	template<typename Matrix, typename Vector> 
	class CholeskyDecomposition<Matrix, Vector, PETSC> : public Factorization<Matrix, Vector, PETSC> {
	public:

		void set_library_type(const DirectSolverLib & TAG)
		{
			this->set_type(TAG, CHOLESKY_DECOMPOSITION_TAG); 
		}

		CholeskyDecomposition()
		{
#ifdef PETSC_HAVE_MUMPS			
			this->set_type(MUMPS_TAG, CHOLESKY_DECOMPOSITION_TAG);
#else
			this->set_type(PETSC_TAG, CHOLESKY_DECOMPOSITION_TAG);
#endif			
		}

	};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	template<typename Matrix, typename Vector> 
	class LUDecomposition<Matrix, Vector, PETSC> : public Factorization<Matrix, Vector, PETSC> {
	public:   
		void set_library_type(const DirectSolverLib & TAG)
		{
			this->set_type(TAG, LU_DECOMPOSITION_TAG); 
		}

		LUDecomposition()
		{
#ifdef PETSC_HAVE_MUMPS			
			this->set_type(MUMPS_TAG, LU_DECOMPOSITION_TAG);
#elif PETSC_HAVE_SUPERLU_DIST
			this->set_type(SUPERLU_DIST_TAG, LU_DECOMPOSITION_TAG);			
#elif PETSC_HAVE_SUPERLU
			this->set_type(SUPERLU_TAG, LU_DECOMPOSITION_TAG);			 // FIXME: runs just with serial matrices 
			                                                      			 // seems that Petsc's also
#else
			this->set_type(PETSC_TAG, LU_DECOMPOSITION_TAG);
#endif			
		}
	};
}





#endif //UTOPIA_PETSC_FACTORIZATIONS_HPP
