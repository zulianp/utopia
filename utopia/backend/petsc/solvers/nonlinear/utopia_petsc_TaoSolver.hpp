#ifndef UTOPIA_PETSC_TAO_SOLVER_HPP
#define UTOPIA_PETSC_TAO_SOLVER_HPP

#include "utopia_BoxConstraints.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_petsc_Types.hpp"
#include <mpi.h>

namespace utopia {
	class TaoSolverWrapper {
	public:
		TaoSolverWrapper();
		~TaoSolverWrapper();
		void destroy();
		bool init(MPI_Comm comm);
		bool set_bounds(const PetscVector &lb, const PetscVector &ub);
		bool solve();

		void set_function(Function<DMatrixd, DVectord> &fun);
		void set_function(Function<DSMatrixd, DVectord> &fun);
		
	private:
		void * data_;
	};


	template<class Matrix, class Vector>
	class TaoSolver final : public NonLinearSolver<Matrix, Vector> {
	public:
		typedef utopia::BoxConstraints<Vector> BoxConstraints;


	private:
		BoxConstraints box_constraints_;
	};
}

#endif //UTOPIA_PETSC_TAO_SOLVER_HPP
