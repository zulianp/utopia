#ifndef UTOPIA_QP_SOLVER_HPP
#define UTOPIA_QP_SOLVER_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_VariableBoundSolverInterface.hpp"

#include <cmath>

namespace utopia 
{

	template<class Matrix, class Vector>
	class QPSolver : 	public IterativeSolver<Matrix, Vector>, 	
						public MatrixFreeLinearSolver<Vector>, 
						public VariableBoundSolverInterface<Vector> 
	{
		public:

			using IterativeSolver<Matrix, Vector>::solve; 
			using MatrixFreeLinearSolver<Vector>::solve; 
	
			QPSolver()
			{

			}		

	        virtual ~QPSolver()
	        {

	        }

			virtual void set_parameters(const Parameters params) override
			{
				IterativeSolver<Matrix, Vector>::set_parameters(params);
			}

			virtual bool solve(const Operator<Vector> &A, const Vector &rhs, Vector &sol) override
			{
				utopia_warning("missing matrix free implementation of give QP-Solver.... \n"); 
				return false; 
			}

	};


}

#endif //UTOPIA_QP_SOLVER_HPP