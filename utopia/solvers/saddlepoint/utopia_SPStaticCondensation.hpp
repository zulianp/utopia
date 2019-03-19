#ifndef UTOPIA_SP_STATIC_CONDENSATION_HPP
#define UTOPIA_SP_STATIC_CONDENSATION_HPP

#include <memory>
#include "utopia_LinearSolver.hpp"
#include "utopia_BiCGStab.hpp"

namespace utopia {

	/**
	@brief Solves the following saddle-point problem:
		| A_m	+ T^T A_s T |  u_m = | f_m + T^Tf_s |
		 u_s  = | T u_m |
	m := master, s := slave
	**/

	template<class Matrix, class Vector>
	class SPStaticCondensation {
	public:

		SPStaticCondensation(const std::shared_ptr<LinearSolver<Matrix, Vector>> &op = std::make_shared<BiCGStab<Matrix, Vector>>())
		: op(op)
		{}

		void linear_solver(const std::shared_ptr<LinearSolver<Matrix, Vector>> &solver)
		{
			op = solver;
		}

		void update(
			const std::shared_ptr<Matrix> &A_m,
			const std::shared_ptr<Matrix> &A_s,
			const std::shared_ptr<Matrix> &T)
		{
			set_up(T);
			update(A_m, A_s);
		}

		void set_up(
			const std::shared_ptr<Matrix> &T)
		{
			this->T = T;
		}

		void update(
			const std::shared_ptr<Matrix> &A_m,
			const std::shared_ptr<Matrix> &A_s)
		{
			this->A_m = A_m;
			this->A_s = A_s;

			S = (*A_m) + transpose(*T) * (*A_s) * (*T);

			op->update(make_ref(S));
		}

		bool apply(
			const Vector &rhs_m,
			const Vector &rhs_s,
			Vector &sol_m,
			Vector &sol_s
		)
		{
			rhs = rhs_m + transpose(*T) * rhs_s;
			bool ok = op->apply(rhs, sol_m);
			sol_s = (*T) * sol_m;
			return ok;
		}

	private:
		Matrix S;
		Vector rhs;
		std::shared_ptr<Matrix> A_m, A_s, T;
		std::shared_ptr<LinearSolver<Matrix, Vector>> op;
	};
}


#endif //UTOPIA_SP_STATIC_CONDENSATION_HPP
