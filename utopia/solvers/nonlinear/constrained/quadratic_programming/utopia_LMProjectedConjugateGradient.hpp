#ifndef UTOPIA_LIMITED_MEMORY_PROJECTED_CONJUGATE_GRADIENT_HPP
#define UTOPIA_LIMITED_MEMORY_PROJECTED_CONJUGATE_GRADIENT_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Smoother.hpp"

#include <cmath>
#include <cassert>

namespace utopia {
	//slow and innefficient implementation just for testing
	template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
	class LMProjectedConjugateGradient : public IterativeSolver<Matrix, Vector> {
	public:
		typedef utopia::BoxConstraints<Vector>  BoxConstraints;
		DEF_UTOPIA_SCALAR(Matrix)

		virtual bool set_box_constraints(const BoxConstraints & box)
		{
			constraints_ = box;
			constraints_.fill_empty_bounds();
			return true;
		}

		virtual void set_parameters(const Parameters params) override
		{
			IterativeSolver<Matrix, Vector>::set_parameters(params);
		}

		bool apply(const Vector &b, Vector &x) override
		{
			return false; 
		}


		bool solve(const Vector &b, Vector &x, const Vector & lb, const Vector & ub) const
		{
			// if(this->verbose()) {
			// 	this->init_solver("LMProjectedConjugateGradient", {" it. ", "|| u - u_old ||"});
			// }

			Vector x_old, x_half, r, uk, wk, zk, pk;

			const SizeType n = local_size(x).get(0); 

			r  = local_zeros(n);
			uk = local_zeros(n);
			wk = local_zeros(n);
			zk = local_zeros(n);
			pk = local_zeros(n);



			// const Matrix &A = *this->get_operator();
			// this->init(local_size(x).get(0)); 
			// const auto &ub = *constraints_.upper_bound();
			// const auto &lb = *constraints_.lower_bound();

			x_old = x;

			Vector Ax; 
			this->apply_H(x, Ax); 

			uk = b - Ax;

			bool converged = false;
			const SizeType check_s_norm_each = 10;			
			pk = -uk;

			int iteration = 0;
			while(!converged) {
				
				// START step
				Vector Ap; 
				this->apply_H(pk, Ap); 

				Scalar alpha = dot(uk, pk)/dot(pk, Ap);
				assert(alpha != 0.);
				if(alpha == 0.) break;
				
				x_half = x_old + alpha * pk;

				x = utopia::min(x_half, ub);
				x = utopia::max(x, lb);

				this->apply_H(x, Ax); 
				uk = b - Ax;

				{
					Write<Vector> w_wk(wk), w_zk(zk);
					Read<Vector> r_uk(uk), r_ub(ub), r_lb(lb), r_p(pk);

					each_read(x, [&](SizeType i, Scalar elem) {
							Scalar val = 0.;
							if (approxeq(elem, ub.get(i)) || approxeq(elem, lb.get(i))) {
								val = std::max(uk.get(i), Scalar(0));
							} else {
								val = uk.get(i);
							}

							if (val == 0) {
								zk.set(i, std::max(pk.get(i), Scalar(0)));
							} else {
								zk.set(i, pk.get(i));
							}

							wk.set(i, val);
	            	});
				}

				this->apply_H(pk, Ap); 
				const Scalar beta = dot(wk, Ap)/dot(pk, Ap);
				pk = wk + beta * zk;

				// END step

				// if(iteration % check_s_norm_each == 0) {
					const Scalar diff = norm2(x_old - x);

					if(this->verbose()) {
					    PrintInfo::print_iter_status({static_cast<Scalar>(iteration), diff});
					}

					//converged = this->check_convergence(iteration, 1, 1, diff);

					if(iteration > 500 || diff < 1e-9)
						converged = true; 
				// }

				++iteration;

				if(converged) break;

				x_old = x;
			}

			return converged;
		}

		// void init(const SizeType & n)
		// {
		// 	r  = local_zeros(n);
		// 	uk = local_zeros(n);
		// 	wk = local_zeros(n);
		// 	zk = local_zeros(n);
		// 	pk = local_zeros(n);
		// }

		virtual void update(const std::shared_ptr<const Matrix> &op) override
		{
		    // IterativeSolver<Matrix, Vector>::update(op);
		    // init(*op);
		}

		LMProjectedConjugateGradient() {}

		LMProjectedConjugateGradient(const LMProjectedConjugateGradient &) = default;

		inline LMProjectedConjugateGradient * clone() const override
		{
			return new LMProjectedConjugateGradient(*this);
		}

		virtual void set_apply_H(std::function< void(const Vector &, Vector &) > fun)
	    {   
	        apply_H = fun; 
	    }

	private:
		BoxConstraints constraints_;
		std::function< void(const Vector &, Vector &) > apply_H; 

		//buffers
		// Vector x_old, x_half, r, uk, wk, zk, pk;
	};
}

#endif //UTOPIA_LIMITED_MEMORY_PROJECTED_CONJUGATE_GRADIENT_HPP
