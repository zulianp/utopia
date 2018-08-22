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
	class ProjectedConjugateGradient : public IterativeSolver<Matrix, Vector>, public Smoother<Matrix, Vector> {
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

		virtual bool smooth(const Vector &b, Vector &x) override
		{
			const Matrix &A = *this->get_operator();

			SizeType it = 0;
			SizeType n_sweeps = this->sweeps();

			while(step(A, b, x) && it++ < n_sweeps) {}

			return it == SizeType(this->sweeps() - 1);
		}

		bool apply(const Vector &b, Vector &x) override
		{
			if(this->verbose())
				this->init_solver("utopia ProjectedConjugateGradient", {" it. ", "|| u - u_old ||"});

			const Matrix &A = *this->get_operator();

			x_old = x;
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

			return true;
		}

		void init(const Matrix &A)
		{

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
		Vector x_old;
	};
}

#endif //UTOPIA_PROJECTED_CONJUGATE_GRADIENT_HPP
