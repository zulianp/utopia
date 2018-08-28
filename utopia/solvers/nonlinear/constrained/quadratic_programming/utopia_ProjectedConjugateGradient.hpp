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

			x_old = x;
			u = b - A * x;

			bool converged = false;
			const SizeType check_s_norm_each = 5;

			int iteration = 0;
			while(!converged) {
				step(A, b, x);

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
			}

			return converged;
		}

		//one step for solving A * x = b : l <= x <= u
		bool step(const Matrix &A, const Vector &b, Vector &x)
		{
			const auto &u = *constraints_.upper_bound();
			const auto &l = *constraints_.lower_bound();

			assert(false && "implement me");


			// Scalar alpha_k_p_1 = dot(u, p)/dot(p, A * p);
			// x_half = x + alpha_k_p_1 * p;
			// x    = max(min(x_half, u), l);
			// u = b - A * x;
			// p = 
			return true;
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
		Vector x_old, x_half, p, u;
	};
}

#endif //UTOPIA_PROJECTED_CONJUGATE_GRADIENT_HPP
