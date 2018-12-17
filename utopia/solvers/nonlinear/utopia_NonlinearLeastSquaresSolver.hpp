#ifndef UTOPIA_NON_LINEAR_LEAST_SQUARES_SOLVER_HPP
#define UTOPIA_NON_LINEAR_LEAST_SQUARES_SOLVER_HPP 

#include "utopia_Function.hpp"
#include "utopia_Parameters.hpp"
#include "utopia_ConvergenceReason.hpp"
#include "utopia_PrintInfo.hpp"
#include "utopia_Monitor.hpp"

namespace utopia 
{
    /**
     * @brief      The base class for all nonlinear solvers. Class provides basic functions used in all nonlinear solvers. 
     *
     * @tparam     Matrix  
     * @tparam     Vector  
     */
	template<class Matrix, class Vector>
	class NonLinearLeastSquaresSolver : public NonLinearSolver<Vector>
	{
	public:
		typedef UTOPIA_SCALAR(Vector)    Scalar;
		typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

		typedef utopia::LinearSolver<Matrix, Vector> Solver;


		NonLinearLeastSquaresSolver(const std::shared_ptr<Solver> &linear_solver)   : 
									linear_solver_(linear_solver)
		{

		}

		virtual ~NonLinearLeastSquaresSolver() {}

		virtual bool solve(LeastSquaresFunction<Matrix, Vector> &fun, Vector &x) = 0;


		void set_linear_solver(const std::shared_ptr<Solver> &linear_solver = std::shared_ptr<Solver>())
		{
			linear_solver_ = linear_solver; 
		}

        virtual void read(Input &in) override
        {
            NonLinearSolver<Vector>::read(in); 

            if(linear_solver_) {
                in.get("linear-solver", *linear_solver_);
            }
        }


	protected:
		inline bool linear_solve(const Matrix &mat, const Vector &rhs, Vector &sol)
		{
			linear_solver_->update(make_ref(mat));
			return linear_solver_->apply(rhs, sol);
		}

	    std::shared_ptr<Solver> linear_solver_;     /*!< Linear solver parameters. */  
	};
}



#endif //UTOPIA_NON_LINEAR_LEAST_SQUARES_SOLVER_HPP
