#ifndef UTOPIA_PETSC_GMRES_HPP
#define UTOPIA_PETSC_GMRES_HPP 


#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_petsc_KSPSolver.hpp"

namespace utopia {
	template<typename Matrix, typename Vector> 
	class GMRES<Matrix, Vector, PETSC> : public KSPSolver<Matrix, Vector, PETSC> {
	public:
		GMRES(const std::string &preconditioner = "jacobi")
		: KSPSolver<Matrix, Vector, PETSC>()
		{
			this->pc_type(preconditioner);
			this->ksp_type("gmres");
		}

		virtual GMRES * clone() const override 
		{
		    return new GMRES(*this);
		}
	};
}

#endif //UTOPIA_PETSC_GMRES_HPP
