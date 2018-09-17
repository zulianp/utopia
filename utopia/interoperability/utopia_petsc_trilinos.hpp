#ifndef UTOPIA_PETSC_TRILINOS_HPP
#define UTOPIA_PETSC_TRILINOS_HPP

#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS
#ifdef WITH_PETSC

#include "utopia_petsc_KSPSolver.hpp"
#include "utopia_petsc_Factorization.hpp"
#include "utopia_CrossBackendLinearSolver.hpp"
#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_petsc_Types.hpp"

namespace utopia {

	template<typename Matrix, typename Vector>
	class KSPSolver<Matrix, Vector, TRILINOS> :
		public CrossBackendLinearSolver<
            Matrix, Vector,
            DSMatrixd, DVectord,
            KSPSolver<DSMatrixd, DVectord, PETSC>
            > {};

	template<typename Matrix, typename Vector>
	class Factorization<Matrix, Vector, TRILINOS> :
		public CrossBackendLinearSolver<
            Matrix, Vector,
            DSMatrixd, DVectord,
            Factorization<DSMatrixd, DVectord, PETSC>
            > {};

}

#endif //WITH_TRILINOS
#endif //WITH_PETSC

#endif //UTOPIA_PETSC_TRILINOS_HPP
