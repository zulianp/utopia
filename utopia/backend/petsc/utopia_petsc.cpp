#include "utopia_petsc.hpp"

//explicit instantiations
namespace utopia {
 // template class Wrapper<PetscSparseMatrix, 2>;

	//petsc linear solvers and smoothers
	template class KSPSolver<DSMatrixd, DVectord>;
	template class ConjugateGradient<DSMatrixd, DVectord>;
	template class GaussSeidel<DSMatrixd, DVectord>;
	

	//petsc non-linear solvers
	template class NonLinearGaussSeidel<DSMatrixd, DVectord>;
	template class SemismoothNewton<DSMatrixd, DVectord, PETSC_EXPERIMENTAL>;

	//
}

