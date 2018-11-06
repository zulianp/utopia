#ifndef UTOPIA_QP_SOLVER_HPP
#define UTOPIA_QP_SOLVER_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_IterativeSolver.hpp"

#include <cmath>

namespace utopia 
{

	template<class Matrix, class Vector>
	class QPSolver : 	public IterativeSolver<Matrix, Vector>, 
						public VariableBoundSolverInterface<Vector> 
	{
		public:
		
			typedef UTOPIA_SCALAR(Vector)                           Scalar;
	        typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;

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



			bool apply(const Vector &b, Vector &x) override = 0; 

			virtual void update(const std::shared_ptr<const Matrix> &op) override = 0; 

			inline QPSolver * clone() const override = 0; 

	};
}

#endif //UTOPIA_QP_SOLVER_HPP