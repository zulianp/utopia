#ifndef UTOPIA_SP_BLOCK_CONJUGATE_GRADIENT_HPP
#define UTOPIA_SP_BLOCK_CONJUGATE_GRADIENT_HPP

#include "utopia_Traits.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_LinearSolverInterfaces.hpp"

#include <iostream>

namespace utopia {

	/**
	@brief Solves the following saddle-point problem:
		| A_m		0		B_t | | sol_m | = | rhs_m |
		| 0			A_s		D_t | | sol_s | = | rhs_s |
		| B 		D 		0   | | lagr  | = | 0     |
	m := master, s := slave
	Warning: rows with boundary conditions in B_t and D_t have to be handled outside.

	*/
	template<class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
	class SPBlockConjugateGradient {
	public:
		using Scalar = UTOPIA_SCALAR(Vector);

		SPBlockConjugateGradient() 
		: op_m(std::make_shared<Factorization<Matrix, Vector>>()),
          op_s(std::make_shared<Factorization<Matrix, Vector>>()),
          verbose_(false), atol_(1e-8), rtol_(1e-8), max_it_(1000)
		{}

		void set_master_solver(const std::shared_ptr< LinearSolver<Matrix, Vector> > &op_m)
		{
			this->op_m = op_m;
		}

		void set_slave_solver(const std::shared_ptr< LinearSolver<Matrix, Vector> > &op_s)
		{
			this->op_s = op_s;
		}

		inline void verbose(const bool verbose)
		{
			verbose_ = verbose;
		}

		inline void atol(const Scalar atol)
		{
			atol_ = atol;
		}

		inline void max_it(const int max_it)
		{
			max_it_ = max_it;
		}

		// inline void rtol(const Scalar rtol)
		// {
		// 	rtol_ = rtol;
		// }

		void update(
			const std::shared_ptr<Matrix> &A_m,
			const std::shared_ptr<Matrix> &A_s,
			const std::shared_ptr<Matrix> &B,
			const std::shared_ptr<Matrix> &D,
			const std::shared_ptr<Matrix> &B_t,
			const std::shared_ptr<Matrix> &D_t)
		{
			set_up(B, D, B_t, D_t);
			update(A_m, A_s);
		}

		void set_up(
			const std::shared_ptr<Matrix> &B,
			const std::shared_ptr<Matrix> &D,
			const std::shared_ptr<Matrix> &B_t,
			const std::shared_ptr<Matrix> &D_t)
		{
			this->B = B;
			this->B_t = B_t;
			this->D = D;
			this->D_t = D_t;
		}

		void update(
			const std::shared_ptr<Matrix> &A_m,
			const std::shared_ptr<Matrix> &A_s)
		{
			op_m->update(A_m);
			op_s->update(A_s);
		}

		bool apply(
			const Vector &rhs_m,
			const Vector &rhs_s,
			Vector &sol_m,
			Vector &sol_s,
			Vector &lagr)
		{
			if(empty(lagr)) {
				lagr = local_zeros(local_size(rhs_s));
			}

			this->residual(
				rhs_m,
				rhs_s,
				lagr,
				r
			);

			bool converged = false;

			double it = 0;
			double rho = 1., rho_1 = 1., beta = 0., alpha = 1., r_norm = 9e9;

			Vector lagr_old = lagr;
			Vector p, q, Ap, r_new, z, z_new;
			while(!converged)
			{
				rho = dot(r, r);

				if(rho == 0.) {
					converged = true;
					break;
				}

				if(it > 0)
				{
					beta = rho/rho_1;
					p = r + beta * p;
				}
				else
				{
					p = r;
				}

				this->apply_op(p, q);
				alpha = rho / dot(p, q);

				lagr += alpha * p;
				r -= alpha * q;

				rho_1 = rho;
				it++;
				r_norm = norm2(lagr - lagr_old);

				converged = r_norm < atol_;
				it++;

				if(!converged && it >= max_it_) {
					break;
				}

				if(verbose_) { std::cout << it << " " << r_norm << std::endl; }

				lagr_old = lagr;
			}

			sol_m = local_zeros(local_size(rhs_m));
			sol_s = local_zeros(local_size(rhs_s));

			op_m->apply(rhs_m - (*B_t) * lagr, sol_m);
			op_s->apply(rhs_s  - (*D_t) * lagr, sol_s);

			return converged;
		}

	private:
		inline void apply_op(const Vector &p, Vector &Ap)
		{
			buff_m = (*B_t) * p;
			buff_s = (*D_t) * p;

			op_m->apply(buff_m, solved_m);
			op_s->apply(buff_s, solved_s);

			Ap = (*B) * solved_m + (*D) * solved_s;
		}

		inline void residual(
			const Vector &rhs_m,
			const Vector &rhs_s,
			const Vector &p,
			Vector &r)
		{
			buff_m = rhs_m - (*B_t) * p;
			buff_s = rhs_s - (*D_t) * p;

			op_m->apply(buff_m, solved_m);
			op_s->apply(buff_s, solved_s);

			r = (*B) * solved_m + (*D) * solved_s;
		}

		std::shared_ptr<Matrix> A_m, A_s;

		std::shared_ptr<Matrix> D, D_t;
		std::shared_ptr<Matrix> B, B_t;

		std::shared_ptr< LinearSolver<Matrix, Vector> > op_m;
		std::shared_ptr< LinearSolver<Matrix, Vector> > op_s;

		Vector buff_m, solved_m, buff_s, solved_s, r;

		bool verbose_;
		Scalar atol_, rtol_;
		int max_it_;
	};

}

#endif //UTOPIA_SP_BLOCK_CONJUGATE_GRADIENT_HPP
