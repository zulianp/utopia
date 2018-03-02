#ifndef UTOPIA_PETSC_GMRES_HPP
#define UTOPIA_PETSC_GMRES_HPP 


#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_petsc_KSPSolver.hpp"

namespace utopia {
	template<typename Matrix, typename Vector> 
	class GMRES<Matrix, Vector, PETSC> : /*public Smoother<Matrix, Vector>,*/ public KSPSolver<Matrix, Vector, PETSC> {
	public:
		GMRES(const Parameters params = Parameters(), const std::string &preconditioner = "jacobi")
		: KSPSolver<Matrix, Vector, PETSC>(params), preconditioner_(preconditioner)
		{
			this->pc_type(preconditioner);
			this->ksp_type("gmres");
		}

		void set_parameters(const Parameters params) override {
			Parameters params_copy = params;
			params_copy.lin_solver_type("gmres");
			params_copy.preconditioner_type(preconditioner_.c_str());
			KSPSolver<Matrix, Vector, PETSC>::set_parameters(params_copy);
		}

		inline void pc_type(const std::string & preconditioner) override
		{
			preconditioner_ = preconditioner;
			KSPSolver<Matrix, Vector, PETSC>::pc_type(preconditioner_);
		}

	 //  // //FIXME use parameters from KSPSolver
		// bool smooth(const Matrix &A, const Vector &rhs, Vector &x) override
		// {
		// 	KSP solver;
		// 	KSPCreate(A.implementation().communicator(), &solver);
		// 	KSPSetFromOptions(solver); 
		// 	KSPSetType(solver, KSPGMRES);
		// 	KSPSetInitialGuessNonzero(solver, PETSC_TRUE);
		// 	KSPSetTolerances(solver, 0., 0., PETSC_DEFAULT, this->sweeps());

		// 	KSPSetOperators(solver, raw_type(A), raw_type(A));

		// 	if(this->verbose()) {
		// 		KSPMonitorSet(
		// 			solver,
		// 			[](KSP, PetscInt iter, PetscReal res, void*) -> PetscErrorCode {
		// 				PrintInfo::print_iter_status({static_cast<PetscReal>(iter), res}); 
		// 				return 0;
		// 			},
		// 			nullptr,
		// 			nullptr);
		// 	}

		// 	KSPSetUp(solver);
		// 	KSPSolve(solver, raw_type(rhs), raw_type(x));
		// 	KSPDestroy(&solver);
		// 	return true;
		// }

	private:
		std::string preconditioner_;
	};
}

#endif //UTOPIA_PETSC_GMRES_HPP
