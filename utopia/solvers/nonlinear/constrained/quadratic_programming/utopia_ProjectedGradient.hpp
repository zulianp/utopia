#ifndef UTOPIA_PROJECTED_GRADIENT_HPP
#define UTOPIA_PROJECTED_GRADIENT_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Smoother.hpp"
#include "utopia_Core.hpp"
#include "utopia_Operations.hpp"

#include <cmath>
#include <cassert>

namespace utopia {
	//slow and innefficient implementation just for testing
	template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
	class ProjectedGradient : public IterativeSolver<Matrix, Vector> {
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
				this->init_solver("utopia ProjectedGradient", {" it. ", "|| u - u_old ||"});

			const Matrix &A = *this->get_operator();
			const auto &upbo = *constraints_.upper_bound();
			const auto &lobo = *constraints_.lower_bound();

			x_old = x;
			u = b - A * x;
			p = u;
			Scalar alpha = 1.;

			bool converged = false;
			const SizeType check_s_norm_each = 20;

			int iteration = 0;
			while(!converged) {
				//perform step
				x_half = x + alpha * p;
				x = utopia::min(upbo, utopia::max(lobo, x_half));
				u = b - A * x;

				{
					Read<Vector>  r_u(u), r_x(x), r_upbo(upbo), r_lobo(lobo);
					Write<Vector> w_(p);

					auto r = range(x);

					for(auto i = r.begin(); i != r.end(); ++i) {
						if(approxeq(x.get(i), upbo.get(i)) || approxeq(x.get(i), lobo.get(i))) {
							p.set(i, 0);
						} else {
							p.set(i, u.get(i));
						}
					}
				}


				if(iteration % check_s_norm_each == 0) {
					const Scalar diff = norm2(x_old - x);

					if(this->verbose()) {
					    PrintInfo::print_iter_status({static_cast<Scalar>(iteration), diff});
					}

					converged = this->check_convergence(iteration, 1, 1, diff);
				}

				++iteration;

				if(converged) break;

				x_old = x;
				alpha = dot(u, p)/dot(p, A * p);
			}

			return converged;
		}


		void init(const Matrix &A)
		{
			auto s = local_size(A);
			p = local_zeros(s.get(0));
		}


		virtual void update(const std::shared_ptr<const Matrix> &op) override
		{
		    IterativeSolver<Matrix, Vector>::update(op);
		    init(*op);
		}

		ProjectedGradient()
		{
		}

		ProjectedGradient(const ProjectedGradient &) = default;

		inline ProjectedGradient * clone() const override
		{
			return new ProjectedGradient(*this);
		}

	private:
		BoxConstraints constraints_;

		//buffers
		Vector x_old, x_half, p, u;
	};
}

#endif //UTOPIA_PROJECTED_GRADIENT_HPP
