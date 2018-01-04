#ifndef UTOPIA_PETSC_FACTORIZATIONS_HPP
#define UTOPIA_PETSC_FACTORIZATIONS_HPP 

#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_petsc_Factorization.hpp"

namespace utopia {

	template<typename Matrix, typename Vector> 
	class CholeskyDecomposition<Matrix, Vector, PETSC> : public LinearSolver<Matrix, Vector> {
	public:
		inline bool apply(const Vector &b, Vector &x) override
		{
		    return strategy_.apply(b, x);
		}  

		inline void update(const std::shared_ptr<const Matrix> &op) override
		{
		    strategy_.update(op);
		}     

		inline void set_parameters(const Parameters params) override
		{
		    LinearSolver<Matrix, Vector>::set_parameters(params);
		    strategy_.set_parameters(params);
		} 

		CholeskyDecomposition()
		{
#ifdef PETSC_HAVE_MUMPS			
			strategy_.set_type(MUMPS_TAG, CHOLESKY_DECOMPOSITION_TAG);
#else
			strategy_.set_type(PETSC_TAG, CHOLESKY_DECOMPOSITION_TAG);
#endif			
		}

	private:
		Factorization<Matrix, Vector> strategy_;
	};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	template<typename Matrix, typename Vector> 
	class LUDecomposition<Matrix, Vector, PETSC> : public LinearSolver<Matrix, Vector> {
	public:
		inline bool apply(const Vector &b, Vector &x) override
		{
		    return strategy_.apply(b, x);
		}  

		inline void update(const std::shared_ptr<const Matrix> &op) override
		{
		    strategy_.update(op);
		}     

		inline void set_parameters(const Parameters params) override
		{
		    LinearSolver<Matrix, Vector>::set_parameters(params);
		    strategy_.set_parameters(params);
		} 

		LUDecomposition()
		{
#ifdef PETSC_HAVE_MUMPS			
			strategy_.set_type(MUMPS_TAG, LU_DECOMPOSITION_TAG);
#elif PETSC_HAVE_SUPERLU_DIST
			strategy_.set_type(SUPERLU_DIST_TAG, LU_DECOMPOSITION_TAG);			
#elif PETSC_HAVE_SUPERLU
			strategy_.set_type(SUPERLU_TAG, LU_DECOMPOSITION_TAG);			 // FIXME: runs just with serial matrices 
			                                                      			 // seems that Petsc's also
#else
			strategy_.set_type(PETSC_TAG, LU_DECOMPOSITION_TAG);
#endif			
		}

	private:
		Factorization<Matrix, Vector> strategy_;
	};
}





#endif //UTOPIA_PETSC_FACTORIZATIONS_HPP
