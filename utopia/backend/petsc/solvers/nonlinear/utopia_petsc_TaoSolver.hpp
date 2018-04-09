#ifndef UTOPIA_PETSC_TAO_SOLVER_HPP
#define UTOPIA_PETSC_TAO_SOLVER_HPP

#include "utopia_BoxConstraints.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_petsc_KSPSolver.hpp"
#include "utopia_petsc_Types.hpp"
#include "utopia_Function.hpp"
#include <mpi.h>
#include <string>

namespace utopia {
	class TaoSolverWrapper {
	public:
		TaoSolverWrapper();
		~TaoSolverWrapper();
		void destroy();
		
		bool init(MPI_Comm comm,
				  const std::string &type,
				  const PetscReal gatol,
				  const PetscReal grtol,
				  const PetscReal gttol,
				  const PetscInt  maxits);

		bool set_bounds(const PetscVector &lb, const PetscVector &ub);
		bool solve(PetscVector &x);

		void set_function(Function<DMatrixd, DVectord> &fun);
		void set_function(Function<DSMatrixd, DVectord> &fun);

		void set_ksp_types(const std::string &ksp, const std::string &pc, const std::string &solver_package);
	private:
		void * data_;
		std::string ksp_type_;
		std::string pc_type_;
		std::string solver_package_;
	};


	template<class Matrix, class Vector>
	class TaoSolver final : public NonLinearSolver<Matrix, Vector> {
	public:
		typedef utopia::BoxConstraints<Vector> BoxConstraints;

		TaoSolver(const std::shared_ptr<LinearSolver<Matrix, Vector>> &linear_solver)
		: NonLinearSolver<Matrix, Vector>(linear_solver)
		{
			this->atol(1e-19);
			this->rtol(1e-12); 
			this->stol(1e-19);
		}

		TaoSolver()
		: NonLinearSolver<Matrix, Vector>(nullptr)
		{
			this->atol(1e-19);
			this->rtol(1e-12); 
			this->stol(1e-19);
		}

		void set_type(const std::string &type)
		{
			type_ = type;
		}

		inline void set_ksp_types(const std::string &ksp, const std::string &pc, const std::string &solver_package)
		{
			impl_.set_ksp_types(ksp, pc, solver_package);
		}

		bool solve(Function<Matrix, Vector> &fun, Vector &x)
		{	
			if(this->linear_solver()) {
				auto ksp_solver = std::dynamic_pointer_cast<KSPSolver<Matrix, Vector>>(this->linear_solver());
				if(ksp_solver) {
					set_ksp_types(ksp_solver->ksp_type(), ksp_solver->pc_type(), ksp_solver->solver_package());
				} else {
					m_utopia_warning_once("> FIXME TaoSolver does not support non KSP based linear solvers");
				}
			}

			impl_.init(
				x.implementation().communicator(),
				type_,
				this->atol(),
				this->rtol(), 
				this->stol(),
				this->max_it()
			);
			
			if(box_constraints_.has_bound()) {
				box_constraints_.fill_empty_bounds();
				impl_.set_bounds(
					box_constraints_.lower_bound()->implementation(),
					box_constraints_.upper_bound()->implementation()
				);
			}

			impl_.set_function(fun);
			return impl_.solve(x.implementation());
		}

		bool set_box_constraints(const BoxConstraints &box_constraints)
		{
			box_constraints_ = box_constraints;
			return true;
		}

	private:
		TaoSolverWrapper impl_;
		BoxConstraints box_constraints_;
		std::string type_;
	};
}

#endif //UTOPIA_PETSC_TAO_SOLVER_HPP
