#ifndef UTOPIA_PROJECTED_CONJUGATE_GRADIENT_HPP
#define UTOPIA_PROJECTED_CONJUGATE_GRADIENT_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Smoother.hpp"

#include <cmath>
#include <cassert>

namespace utopia {
	//slow and innefficient implementation just for testing
	template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
	class ProjectedConjugateGradient : public IterativeSolver<Matrix, Vector> {
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
			if(this->verbose())
				this->init_solver("utopia ProjectedConjugateGradient", {" it. ", "|| u - u_old ||"});

			const Matrix &A = *this->get_operator();
			const auto &ub = *constraints_.upper_bound();
			const auto &lb = *constraints_.lower_bound();

			x_old = x;
			u = b - A * x;

			bool converged = false;
			const SizeType check_s_norm_each = 10;			
			p = -u;

			int iteration = 0;
			while(!converged) {
				// START step
				Scalar alpha = dot(u, p)/dot(p, A * p);
				x_half = x_old + alpha * p;
				// std::cout << x_half.size() << std::endl;
				x = utopia::max(utopia::min(x_half, ub), lb);
				uk = b - A * x;

				{
					Write<Vector> w_wk(wk), w_zk(zk);
					Read<Vector> r_uk(uk), r_lb(lb), r_p(p);

					each_read(x, [&](SizeType i, Scalar elem) {
							Scalar val = 0.;
							if (approxeq(elem, ub.get(i)) || approxeq(elem, lb.get(i))) {
								val = std::max(uk.get(i), Scalar(0));
							} else {
								val = uk.get(i);
							}

							if (val == 0) {
								zk.set(i, std::max(p.get(i), Scalar(0)));
							} else {
								zk.set(i, p.get(i));
							}

							wk.set(i, val);
	            });
				}


				Scalar beta = dot(wk, A * p)/dot(p, A * p);
				pk = wk + beta * zk;

				// END step

				if(iteration % check_s_norm_each == 0) {
					const Scalar diff = norm2(x_old - x);

					if(this->verbose()) {
					    PrintInfo::print_iter_status({static_cast<Scalar>(iteration), diff});
					}

					converged = this->check_convergence(iteration, 1, 1, diff);
				}

				++iteration;

				if(converged) break;

				p = pk;
        		u = uk;
				x_old = x;
			}

			return converged;
		}

		//one step for solving A * x = b : l <= x <= u
		// bool step(const Matrix &A, const Vector &b, Vector &x)
		// {
		// 	const auto &u = *constraints_.upper_bound();
		// 	const auto &l = *constraints_.lower_bound();
		//
		// 	assert(false && "implement me");
		//
		//
		// 	// Scalar alpha_k_p_1 = dot(r, p)/dot(p, A * p);
		// 	// x_half = x + alpha_k_p_1 * p;
		// 	// x    = max(min(x_half, u), l);
		// 	// u = b - A * x;
		// 	// p =
		// 	return true;
		// }

		void init(const Matrix &A)
		{
			auto n = local_size(A).get(0);
			r  = local_zeros(n);
			uk = local_zeros(n);
			wk = local_zeros(n);
			zk = local_zeros(n);
			pk = local_zeros(n);
			p  = local_zeros(n);
		}


		virtual void update(const std::shared_ptr<const Matrix> &op) override
		{
		    IterativeSolver<Matrix, Vector>::update(op);
		    init(*op);
		}

		ProjectedConjugateGradient()
		{
		}

		ProjectedConjugateGradient(const ProjectedConjugateGradient &) = default;

		inline ProjectedConjugateGradient * clone() const override
		{
			return new ProjectedConjugateGradient(*this);
		}

	private:
		BoxConstraints constraints_;

		//buffers
		Vector x_old, x_half, u, r, uk, wk, zk, pk, p;
	};
}

#endif //UTOPIA_PROJECTED_CONJUGATE_GRADIENT_HPP
