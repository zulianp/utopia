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
						public virtual VariableBoundSolverInterface<Vector> 
	{
		public:	
			QPSolver()
			{

			}		

	        virtual ~QPSolver()
	        {

	        }

	        virtual QPSolver * clone() const override = 0; 
	};


	template<class Vector>
	class MatrixFreeQPSolver : 	public MatrixFreeLinearSolver<Vector>, 
								public virtual VariableBoundSolverInterface<Vector> 
	{
		public:	
			MatrixFreeQPSolver()
			{

			}		

	        virtual ~MatrixFreeQPSolver()
	        {

	        }


	        virtual MatrixFreeQPSolver * clone() const override = 0; 

	};


}

#endif //UTOPIA_QP_SOLVER_HPP