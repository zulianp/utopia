#ifndef UTOPIA_PETSC_TAO_SOLVER_HPP
#define UTOPIA_PETSC_TAO_SOLVER_HPP

#include "utopia_BoxConstraints.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_petsc_Types.hpp"
#include "utopia_petsc_Function.hpp"
#include <mpi.h>

namespace utopia {
	class TaoSolverWrapper {
	public:
		TaoSolverWrapper();
		~TaoSolverWrapper();
		void destroy();
		bool init(MPI_Comm comm);
		bool set_bounds(const PetscVector &lb, const PetscVector &ub);
		bool solve(PetscVector &x);

		void set_function(Function<DMatrixd, DVectord> &fun);
		void set_function(Function<DSMatrixd, DVectord> &fun);
		
	private:
		void * data_;
	};


	template<class Matrix, class Vector>
	class TaoSolver final : public NonLinearSolver<Matrix, Vector> {
	public:
		typedef utopia::BoxConstraints<Vector> BoxConstraints;

		TaoSolver(const std::shared_ptr<LinearSolver<Matrix, Vector>> &linear_solver)
		: NonLinearSolver<Matrix, Vector>(linear_solver)
		{}

		virtual bool solve(Function<Matrix, Vector> &fun, Vector &x)
		{
			impl_.init(x.implementation().communicator());
			
			if(box_constraints_.has_bound()) {
				impl_.set_bounds(
					box_constraints_.lower_bound()->implementation(),
					box_constraints_.upper_bound()->implementation()
				);
			}

			impl_.set_function(fun);
			return impl_.solve(x.implementation());
		}

	private:
		TaoSolverWrapper impl_;
		BoxConstraints box_constraints_;
	};
}

#endif //UTOPIA_PETSC_TAO_SOLVER_HPP
