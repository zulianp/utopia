#ifndef UTOPIA_PETSC_FACTORIZATIONS_HPP
#define UTOPIA_PETSC_FACTORIZATIONS_HPP 

#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_petsc_Factorization.hpp"
#include "utopia_SolverType.hpp"

namespace utopia {

	template<typename Matrix, typename Vector> 
	class CholeskyDecomposition<Matrix, Vector, PETSC> final: public Factorization<Matrix, Vector, PETSC> {
	public:

		void set_library_type(const SolverPackage &package)
		{
			this->set_type(package, Solver::cholesky_decomposition()); 
		}

		CholeskyDecomposition()
		{
#ifdef PETSC_HAVE_MUMPS			
			this->set_type(Solver::mumps(), Solver::cholesky_decomposition());
#else
			this->set_type(Solver::petsc(), Solver::cholesky_decomposition());
#endif			
		}

	};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	template<typename Matrix, typename Vector> 
	class LUDecomposition<Matrix, Vector, PETSC> final: public Factorization<Matrix, Vector, PETSC> {
	public:   
		void set_library_type(const SolverPackage &package)
		{
			this->set_type(package, Solver::lu_decomposition()); 
		}

		LUDecomposition()
		{
#ifdef PETSC_HAVE_MUMPS			
			this->set_type(Solver::mumps(), Solver::lu_decomposition());
#elif PETSC_HAVE_SUPERLU_DIST
			this->set_type(Solver::superlu_dist(), Solver::lu_decomposition());			
#elif PETSC_HAVE_SUPERLU
			this->set_type(Solver::superlu(), Solver::lu_decomposition());			 // FIXME: runs just with serial matrices 
			                                                      			 // seems that Petsc's also
#else
			this->set_type(Solver::petsc(), Solver::lu_decomposition());
#endif			
		}
	};
}





#endif //UTOPIA_PETSC_FACTORIZATIONS_HPP
