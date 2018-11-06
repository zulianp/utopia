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
						public VariableBoundSolverInterface<Vector> 
	{
		public:
	
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

	};
}

#endif //UTOPIA_QP_SOLVER_HPP