/*
 * @Author: alenakopanicakova
 * @Date:   2016-06-06
 * @Last Modified by:   alenakopanicakova
 * @Last Modified time: 2016-11-06
 */
#ifndef UTOPIA_CONJUGATE_GRAD_H
#define UTOPIA_CONJUGATE_GRAD_H

#include "utopia_IterativeSolver.hpp"
#include "utopia_MatrixFreeLinearSolver.hpp"
#include "utopia_Parameters.hpp"
#include "utopia_Preconditioner.hpp"
#include "utopia_Smoother.hpp"

#include <memory>

namespace utopia 
{
	

	//FIXME also use the PreconditionedSolver interface properly
	/**
	 * @brief      Conjugate Gradient solver. Works with all utopia tensor types.
	 * @tparam     Matrix
	 * @tparam     Vector
	 */
	template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
	class ConjugateGradient : public IterativeSolver<Matrix, Vector>, public Smoother<Matrix, Vector>, public MatrixFreeLinearSolver<Vector>
	{
		typedef UTOPIA_SCALAR(Vector) 	 Scalar;
		typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
		typedef utopia::LinearSolver<Matrix, Vector> Solver;
		typedef utopia::Preconditioner<Vector> Preconditioner;

	public:

		using IterativeSolver<Matrix, Vector>::solve;
		
		ConjugateGradient(const Parameters params = Parameters())
		: reset_initial_guess_(false)
		{
			set_parameters(params);
		}
		
		void reset_initial_guess(const bool val)
		{
			reset_initial_guess_ = val;
		}

		/**
		 * @brief      Sets the parameters.
		 *
		 * @param[in]  params  The parameters
		 */
		void set_parameters(const Parameters params) override
		{
			IterativeSolver<Matrix, Vector>::set_parameters(params);
		}
		
		/**
		 * @brief      Solution routine for CG.
		 *
		 * @param[in]  b     The right hand side.
		 * @param      x     The initial guess/solution.
		 *
		 * @return true if the linear system has been solved up to required tollerance. False otherwise
		 */
		bool apply(const Vector &b, Vector &x) override
		{
			auto A_ptr = utopia::op(this->get_operator());
			return solve(*A_ptr, b, x);
		}

		
		bool solve(const Operator<Vector> &A, const Vector &b, Vector &x) override
		{
			init(local_size(b).get(0));

			if(precond_) {
				return preconditioned_solve(A, b, x);
			} else {
				return unpreconditioned_solve(A, b, x);
			}
		}

		/**
		 * @brief      Sets the preconditioner.
		 *
		 * @param[in]  precond  The preconditioner.
		 */
		void set_preconditioner(const std::shared_ptr<Preconditioner> &precond)
		{
			precond_ = precond;
		}
		
		/*! @brief if overriden the subclass has to also call this one first
		 */
		virtual void update(const std::shared_ptr<const Matrix> &op) override
		{
			IterativeSolver<Matrix, Vector>::update(op);

			init(local_size(*op).get(0));
			
			if(precond_) {
				auto ls_ptr = dynamic_cast<LinearSolver<Matrix, Vector> *>(precond_.get());
				if(ls_ptr) {
					ls_ptr->update(op);
				}
			}
		}

		bool smooth(const Vector &rhs, Vector &x) override
		{
			SizeType temp = this->max_it();
			this->max_it(this->sweeps());
			auto A_ptr = utopia::op(this->get_operator());
			unpreconditioned_solve(*A_ptr, rhs, x);
			this->max_it(temp);
			return true;
		}
		
		ConjugateGradient * clone() const override
		{
			return new ConjugateGradient(*this);
		}
		
	private:
		bool unpreconditioned_solve(const Operator<Vector> &A, const Vector &b, Vector &x)
		{
			SizeType it = 0;
			Scalar rho = 1., rho_1 = 1., beta = 0., alpha = 1., r_norm = 9e9;
			
			assert(!empty(b));
			
			if(empty(x) || size(x).get(0) != size(b).get(0)) {
				x = local_zeros(local_size(b));
				r = b;
			} else {
				assert(local_size(x).get(0) == local_size(b).get(0));
				if(reset_initial_guess_) {
					x.set(0.);
				}
				// r = b - A * x;
				A.apply(x, r);
				r = b - r;
			}
			
			this->init_solver("Utopia Conjugate Gradient", {"it. ", "||r||" });
			bool converged = false;

			SizeType check_norm_each = 1;
			
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
				
				// q = A * p;
				A.apply(p, q);
				alpha = rho / dot(p, q);
				
				x += alpha * p;
				r -= alpha * q;
				
				rho_1 = rho;

				if((it % check_norm_each) == 0) {
					// r = 
					// A.apply(x, r);
					// r = b - r;
					
					r_norm = norm2(r);
					
					if(this->verbose()) {
						PrintInfo::print_iter_status({it, r_norm });
					}
					
					converged = this->check_convergence(it, r_norm, 1, 1);
				}

				it++;
			}
			
			return converged;
		}
		
		bool preconditioned_solve(const Operator<Vector> &A, const Vector &b, Vector &x)
		{
			SizeType it = 0;
			Scalar beta = 0., alpha = 1., r_norm = 9e9;
			
			z     = local_zeros(local_size(b));
			z_new = local_zeros(local_size(b));
			
			if(empty(x) || size(x).get(0) != size(b).get(0)) {
				x = local_zeros(local_size(b));
				r = b;
			} else {
				assert(local_size(x).get(0) == local_size(b).get(0));
				
				if(reset_initial_guess_) {
					x.set(0.);
				}

				A.apply(x, r);
				r = b - r;
			}
			
			precond_->apply(r, z);
			p = z;
			
			this->init_solver("Utopia Conjugate Gradient", {"it. ", "||r||" });
			bool stop = false;
			
			while(!stop)
			{
				// Ap = A*p;
				A.apply(p, Ap);
				alpha = dot(r, z)/dot(p, Ap);

				if(std::isinf(alpha) || std::isnan(alpha)) {
					stop = this->check_convergence(it, r_norm, 1, 1);
					break;
				}
				
				x += alpha * p;
				r_new = r - alpha * Ap;
				
				r_norm = norm2(r_new);
				
				if(r_norm < this->atol()) {
					if(this->verbose()) {
						PrintInfo::print_iter_status({it, r_norm});
					}
					
					stop = this->check_convergence(it, r_norm, 1, 1);
					break;
				}
				
				precond_->apply(r_new, z_new);
				beta = dot(z_new, r_new)/dot(z, r);
				
				p = z_new + beta * p;
				r = r_new;
				z = z_new;
				
				if(this->verbose()) {
					PrintInfo::print_iter_status({it, r_norm});
				}
				
				stop = this->check_convergence(it, r_norm, 1, 1);
				it++;
			}
			
			if(r_norm <= this->atol()) {
				//FIXME sometimes this fails for some reason
				// assert(check_solution(A, x, b));
				return true;
			} else {
				return false;
			}
		}

		bool check_solution(const Operator<Vector> &A, const Vector &x, const Vector &b) const
		{
			Vector r;
			A.apply(x, r);
			r -= b;

			const Scalar r_norm = norm2(r);
			
			if(r_norm > 100 * this->atol()) {
				// write("A.m", *this->get_operator());
				// disp(*this->get_operator());
				assert(r_norm <= this->atol());
				return false;
			}

			return true;
		}

		void init(const SizeType &ls)
		{
			auto zero_expr = local_zeros(ls);

			//resets all buffers in case the size has changed
			if(!empty(r)) {
				r = zero_expr;
			}

			if(!empty(p)) {
				p = zero_expr;
			}

			if(!empty(q)) {
				q = zero_expr;
			}

			if(!empty(Ap)) {
				Ap = zero_expr;
			}

			if(!empty(r_new)) {
				r_new = zero_expr;
			}

			if(!empty(z)) {
				z = zero_expr;
			}

			if(!empty(z_new)) {
				z_new = zero_expr;
			}
		}
		
		std::shared_ptr<Preconditioner> precond_;
		Vector r, p, q, Ap, r_new, z, z_new;
		bool reset_initial_guess_;
	};
}


#endif //UTOPIA_CONJUGATE_GRAD_H

