#ifndef UTOPIA_PROJECTED_CONJUGATE_GRADIENT_HPP
#define UTOPIA_PROJECTED_CONJUGATE_GRADIENT_HPP

#include "utopia_ForwardDeclarations.hpp"


#include <cmath>
#include <cassert>

namespace utopia 
{
	//slow and innefficient implementation just for testing
	template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
	class ProjectedConjugateGradient : public QPSolver<Matrix, Vector> {
	public:

		DEF_UTOPIA_SCALAR(Matrix)

		ProjectedConjugateGradient() {}

		ProjectedConjugateGradient(const ProjectedConjugateGradient &) = default;

		inline ProjectedConjugateGradient * clone() const override
		{
			auto ptr = new ProjectedConjugateGradient(*this);
			ptr->set_box_constraints(this->get_box_constraints());
			return ptr;
		}

		virtual bool apply(const Vector &b, Vector &x) override
		{
			if(this->verbose()) {
				this->init_solver("ProjectedConjugateGradient", {" it. ", "|| u - u_old ||"});
			}

			const Matrix &A = *this->get_operator();

			// ideally, we have two separate implementations, or cases
			this->fill_empty_bounds(); 
			
			const auto &ub = this->get_upper_bound();
			const auto &lb = this->get_lower_bound();			

			x_old = x;
			uk = b - A * x;

			bool converged = false;
			const SizeType check_s_norm_each = 10;			
			pk = -uk;

			int iteration = 0;
			while(!converged) {
				
				// START step

				Scalar alpha = dot(uk, pk)/dot(pk, A * pk);
				assert(alpha != 0.);
				if(alpha == 0. || std::isinf(alpha) || std::isnan(alpha)) break;
				
				x_half = x_old + alpha * pk;

				x = utopia::min(x_half, ub);
				x = utopia::max(x, lb);

				uk = b - A * x;

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

				const Scalar beta = dot(wk, A * pk)/dot(pk, A * pk);
				pk = wk + beta * zk;

				// END step

				if(iteration % check_s_norm_each == 0 || std::isinf(beta) || std::isnan(beta)) {
					const Scalar diff = norm2(x_old - x);

					if(this->verbose()) {
					    PrintInfo::print_iter_status({static_cast<Scalar>(iteration), diff});
					}

					converged = this->check_convergence(iteration, 1, 1, diff);
				}

				++iteration;

				if(converged) break;

				x_old = x;
			}

			return converged;
		}

		void init(const Matrix &A)
		{
			auto n = local_size(A).get(0);
			r  = local_zeros(n);
			uk = local_zeros(n);
			wk = local_zeros(n);
			zk = local_zeros(n);
			pk = local_zeros(n);
		}

		virtual void update(const std::shared_ptr<const Matrix> &op) override
		{
		    QPSolver<Matrix, Vector>::update(op);
		    init(*op);
		}


	private:
		//buffers
		Vector x_old, x_half, r, uk, wk, zk, pk;
	};
}

#endif //UTOPIA_PROJECTED_CONJUGATE_GRADIENT_HPP
