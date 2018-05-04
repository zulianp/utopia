#ifndef UTOPIA_PETSC_GMRES_HPP
#define UTOPIA_PETSC_GMRES_HPP 


#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_petsc_KSPSolver.hpp"

namespace utopia {
	template<typename Matrix, typename Vector> 
	class GMRES<Matrix, Vector, PETSC> : public KSPSolver<Matrix, Vector, PETSC> {
	public:
		GMRES(const Parameters params = Parameters(), const std::string &preconditioner = "jacobi")
		: KSPSolver<Matrix, Vector, PETSC>(params)
		{
			this->pc_type(preconditioner);
			this->ksp_type("gmres");
		}

		void set_parameters(const Parameters params) override {
			Parameters params_copy = params;
			params_copy.lin_solver_type("gmres");
			params_copy.preconditioner_type(this->pc_type().c_str());
			KSPSolver<Matrix, Vector, PETSC>::set_parameters(params_copy);
		}

		virtual GMRES * clone() const override 
		{
		    return new GMRES(*this);
		}
	};
}

#endif //UTOPIA_PETSC_GMRES_HPP
