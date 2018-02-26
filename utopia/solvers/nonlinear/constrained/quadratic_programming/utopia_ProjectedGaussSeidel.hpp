#ifndef UTOPIA_PROJECTED_GAUSS_SEIDEL_HPP
#define UTOPIA_PROJECTED_GAUSS_SEIDEL_HPP 

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Smoother.hpp"


namespace utopia {
	//slow and innefficient implementation just for testing
	template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
	class ProjectedGaussSeidel : public IterativeSolver<Matrix, Vector>, public Smoother<Matrix, Vector> {
	public:
		typedef utopia::BoxConstraints<Vector>  BoxConstraints;
		DEF_UTOPIA_SCALAR(Matrix)

		virtual bool set_box_constraints(const BoxConstraints & box)
		{
			constraints_ = box;
			return true;
		}

		virtual void set_parameters(const Parameters params) override
		{
			IterativeSolver<Matrix, Vector>::set_parameters(params);
		}

		virtual bool smooth(const Matrix &A, const Vector &b, Vector &x) override
		{
			init(A);
			std::size_t it = 0;
			if(constraints_.has_bound()) {
				while(step(A, b, x) && it++ < this->sweeps()) {}
			} else {
				while(unconstrained_step(A, b, x) && it++ < this->sweeps()) {}
			}
			return it == this->sweeps() - 1;
		}

		bool apply(const Vector &b, Vector &x) override
		{ 
			if(this->verbose())
				this->init_solver("utopia ProjectedGaussSeidel", {" it. ", "|| u - u_old ||"});

			const Matrix &A = *this->get_operator();

			//TODO generic version
			assert(!constraints_.has_lower_bound() );
			init(A);

			x_old = x;
			bool converged = false;

			int iteration = 0;
			while(!converged) {
				if(constraints_.has_bound()) {
					step(A, b, x);
				} else {
					unconstrained_step(A, b, x);
				}

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

		void non_linear_jacobi_step(const Matrix &A, const Vector &b, Vector &x)
		{
			r = b - A * x;
			x = min(x + e_mul(d_inv, r), *constraints_.upper_bound());
		}

		bool unconstrained_step(const Matrix &A, const Vector &b, Vector &x)
		{
			r = b - A * x;
			c *= 0.;
			
			Range rr = row_range(A);
			{
				ReadAndWrite<Vector> rw_c(c);
				Read<Vector> r_r(r), r_d_inv(d_inv);
				Read<Matrix> r_A(A);

				for(SizeType il = 0; il < this->n_local_sweeps(); ++il) {
					for(auto i = rr.begin(); i != rr.end(); ++i) {
						RowView<const Matrix> row_view(A, i);

						auto s = r.get(i);

						for(auto index = 0; index < row_view.n_values(); ++index) {
							const auto j    = row_view.col(index);
							const auto a_ij = row_view.get(index);

							if(rr.inside(j) && i != j) {
								s -= a_ij * c.get(j);
							}
						}

						//update correction
						c.set(i, d_inv.get(i) * s );
					}

					if(use_symmetric_sweep_) {
						for(auto i = rr.end()-1; i >= rr.begin(); --i) {
							RowView<const Matrix> row_view(A, i);

							auto s = r.get(i);

							for(auto index = 0; index < row_view.n_values(); ++index) {
								const auto j    = row_view.col(index);
								const auto a_ij = row_view.get(index);

								if(rr.inside(j) && i != j) {
									s -= a_ij * c.get(j);
								}
							}

							//update correction
							c.set(i, d_inv.get(i) * s);
						}
					}
				}
			}

			Scalar alpha = 1.;

			if(use_line_search_) {
				const Scalar rho = dot(c, r);
				const Scalar denom = dot(A * c, c);

				if(denom == 0.) return true;

				alpha = rho/denom;

				assert(alpha > 0);

				if(alpha <= 0) {
					std::cerr << "[Warning] negative alpha" << std::endl;
					alpha = 1.;
					c = r;
				}
			}

			x += alpha * c;
			return true;
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
				Read<Vector> r_r(r), r_d_inv(d_inv), r_g(g);
				Read<Matrix> r_A(A);

				for(SizeType il = 0; il < this->n_local_sweeps(); ++il) {

					for(auto i = rr.begin(); i != rr.end(); ++i) {
						RowView<const Matrix> row_view(A, i);

						auto s = r.get(i);

						for(auto index = 0; index < row_view.n_values(); ++index) {
							const auto j    = row_view.col(index);
							const auto a_ij = row_view.get(index);

							if(rr.inside(j) && i != j) {
								s -= a_ij * c.get(j);
							}
						}

						//update correction
						c.set(i, std::min( d_inv.get(i) * s, g.get(i)) );
					}
				
					if(use_symmetric_sweep_) {
						for(auto i = rr.end()-1; i >= rr.begin(); --i) {
							RowView<const Matrix> row_view(A, i);

							auto s = r.get(i);

							for(auto index = 0; index < row_view.n_values(); ++index) {
								const auto j    = row_view.col(index);
								const auto a_ij = row_view.get(index);

								if(rr.inside(j) && i != j) {
									s -= a_ij * c.get(j);
								}
							}

							//update correction
							c.set(i, std::min( d_inv.get(i) * s, g.get(i)) );
						}
					}
				}
			}

			Scalar alpha = 1.;
			
			if(use_line_search_) {
				inactive_set_ *= 0.;

				{
					Read<Vector> r_c(c), r_g(g);
					Write<Vector> w_a(inactive_set_);

					for(auto i = rr.begin(); i != rr.end(); ++i) {
						if(c.get(i) < g.get(i)) {
							inactive_set_.set(i, 1.);
						}
					}
				}

				is_c_ = e_mul(c, inactive_set_);

				const Scalar rho = dot(is_c_, r);
				alpha = rho/dot(A * is_c_, is_c_);

				assert(alpha > 0);
				
				if(alpha <= 0) {
					std::cerr << "[Warning] negative alpha" << std::endl;
					alpha = 1.;
					descent_dir = r;
				} else {
					descent_dir = alpha * c;
				}
			}

			x += utopia::min(descent_dir, g);
			return true;
		}

		void init(const Matrix &A)
		{
			d = diag(A);
			d_inv = 1./d;
			c = local_zeros(local_size(A).get(0));

			if(use_line_search_) {
				inactive_set_ = local_zeros(local_size(c));
			}
		}


		virtual void update(const std::shared_ptr<const Matrix> &op) override
		{
		    IterativeSolver<Matrix, Vector>::update(op);
		    init(*op);
		}

		ProjectedGaussSeidel()
		: use_line_search_(true), use_symmetric_sweep_(true), n_local_sweeps_(3)
		{}

		void set_use_line_search(const bool val) 
		{
			use_line_search_ = val;
		}

		inline SizeType n_local_sweeps() const
		{
			return n_local_sweeps_;
		}

		inline void set_n_local_sweeps(const SizeType n_local_sweeps)
		{
			n_local_sweeps_ = n_local_sweeps;
		}

		inline void set_use_symmetric_sweep(const bool use_symmetric_sweep)
		{
			use_symmetric_sweep_ = use_symmetric_sweep;
		}


		
	private:
		bool use_line_search_;
		bool use_symmetric_sweep_;
		SizeType n_local_sweeps_;

		BoxConstraints constraints_;	

		Vector r, d, g, c, d_inv, x_old, descent_dir;
		Vector inactive_set_;
		Vector is_c_;
	};
}

#endif //UTOPIA_PROJECTED_GAUSS_SEIDEL_HPP
