#ifndef UTOPIA_PROJECTED_GAUSS_SEIDEL_HPP
#define UTOPIA_PROJECTED_GAUSS_SEIDEL_HPP 

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_IterativeSolver.hpp"

namespace utopia {
	//slow and innefficient implementation just for testing
	template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
	class ProjectedGaussSeidel : public IterativeSolver<Matrix, Vector> {
	public:
		typedef utopia::BoxConstraints<Vector>  BoxConstraints;
		DEF_UTOPIA_SCALAR(Matrix);

		virtual bool set_box_constraints(const BoxConstraints & box)
		{
			constraints_ = box;
			return true;
		}

		virtual void set_parameters(const Parameters params) override
		{
			IterativeSolver<Matrix, Vector>::set_parameters(params);
		}

		bool solve(const Matrix &A, const Vector &b, Vector &x) override
		{ 
			if(this->verbose())
				this->init_solver("utopia ProjectedGaussSeidel", {" it. ", "|| u - u_old ||"});

			//TODO generic version
			assert( constraints_.has_upper_bound() && !constraints_.has_lower_bound() );
			init(A);

			x_old = x;
			bool converged = false;

			int iteration = 0;
			while(!converged) {
				step(A, b, x);
				const Scalar diff = norm2(x_old - x);

				if(this->verbose())
				    PrintInfo::print_iter_status({static_cast<Scalar>(iteration), diff}); 

				converged = this->check_convergence(iteration, 1, 1, diff);


				++iteration;
				if(converged) break;
				x_old = x;
			}

			return converged;
		}

		bool step(const Matrix &A, const Vector &b, Vector &x)
		{
			r = b - A * x;
			//localize gap function for correction
			g = *constraints_.upper_bound() - x;
			c *= 0.;
			
			Range rr = row_range(A);
			{
				ReadAndWrite<Vector> rw_c(c);
				Read<Vector> r_r(r), r_d_inv(d_inv);

				for(auto i = rr.begin(); i != rr.end(); ++i) {
					RowView<const Matrix> row_view(A, i);

					auto s = r.get(i);

					for(auto index = 0; index < row_view.n_values(); ++index) {
						const auto j   = row_view.get_col_at(index);
						const auto val = row_view.get_value_at(index);

						if(rr.inside(j) && i != j) {
							s -= val * c.get(j);
						}
					}

					//update correction
					c.set(i, std::min( d_inv.get(i) * s, g.get(i)) );
				}
			}

			x += c;
			return true;
		}

		void init(const Matrix &A)
		{
			d = diag(A);
			d_inv = 1./d;
			c = local_zeros(local_size(A).get(0));
		}
		
	private:
		BoxConstraints constraints_;	
		Vector r, d, g, c, d_inv, x_old;
	};
}

#endif //UTOPIA_PROJECTED_GAUSS_SEIDEL_HPP
