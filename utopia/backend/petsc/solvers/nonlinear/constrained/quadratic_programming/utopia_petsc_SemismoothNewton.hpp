#ifndef UTOPIA_PETSC_SEMI_SMOOTH_NEWTON_HPP
#define UTOPIA_PETSC_SEMI_SMOOTH_NEWTON_HPP

#include "utopia_SemismoothNewton.hpp"
#include "utopia_petsc.hpp"
#include "utopia_QPSolver.hpp"
#include <memory>

namespace utopia {

	template<class Matrix, class Vector>
	class SemismoothNewton<Matrix, Vector, PETSC_EXPERIMENTAL> final : public QPSolver<Matrix, Vector> {
		typedef UTOPIA_SCALAR(Vector)    			 Scalar;
		typedef UTOPIA_SIZE_TYPE(Vector) 			 SizeType;
		typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;
		typedef utopia::QPSolver<Matrix, Vector> 	 Super;
		typedef utopia::BoxConstraints<Vector>       BoxConstraints;

	public:

		SemismoothNewton(const std::shared_ptr<LinearSolver> &linear_solver = std::make_shared<Factorization<Matrix, Vector>>());
		~SemismoothNewton();
		SemismoothNewton * clone() const override;
		bool apply(const Vector &rhs, Vector &sol) override;
		void read(Input &in) override;

	private:
		class Impl;
		std::unique_ptr<Impl> impl_;

	};

}

#endif //UTOPIA_PETSC_SEMI_SMOOTH_NEWTON_HPP
