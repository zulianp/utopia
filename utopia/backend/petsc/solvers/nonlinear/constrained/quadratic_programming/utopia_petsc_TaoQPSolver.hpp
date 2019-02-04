#ifndef UTOPIA_PETSC_TAO_QP_SOLVEWR_HPP
#define UTOPIA_PETSC_TAO_QP_SOLVEWR_HPP

#include "utopia_QPSolver.hpp"
#include "utopia_LinearSolverInterfaces.hpp"
#include <memory>

namespace utopia {

	template<class Matrix, class Vector>
	class TaoQPSolver final : public QPSolver<Matrix, Vector> {
		typedef UTOPIA_SCALAR(Vector)    			 Scalar;
		typedef UTOPIA_SIZE_TYPE(Vector) 			 SizeType;
		typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;
		typedef utopia::QPSolver<Matrix, Vector> 	 Super;
		typedef utopia::BoxConstraints<Vector>       BoxConstraints;

	public:

		TaoQPSolver(const std::shared_ptr<LinearSolver> &linear_solver = std::make_shared<Factorization<Matrix, Vector>>());
		~TaoQPSolver();
		TaoQPSolver * clone() const override;
		bool apply(const Vector &rhs, Vector &sol) override;
		void pc_type(const std::string & pc_type);
		void tao_type(const std::string &type);
	private:
		class Impl;
		std::unique_ptr<Impl> impl_;

	};

}

#endif //UTOPIA_PETSC_TAO_QP_SOLVEWR_HPP
