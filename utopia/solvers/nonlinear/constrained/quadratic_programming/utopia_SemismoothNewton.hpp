#ifndef UTOPIA_SOLVER_SEMISMOOTH_NEWTON_HPP
#define UTOPIA_SOLVER_SEMISMOOTH_NEWTON_HPP

#include "utopia_Wrapper.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Core.hpp"
#include <vector>
#include <memory>

namespace utopia {
	
	template<class Matrix, class Vector>
	class SemismoothNewton : public IterativeSolver<Matrix, Vector> {
		typedef UTOPIA_SCALAR(Vector)    Scalar;
		typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
		typedef utopia::LinearSolver<Matrix, Vector> Solver;
		typedef utopia::BoxConstraints<Vector>      BoxConstraints;
		
	public:
		
		SemismoothNewton(const std::shared_ptr <Solver> &linear_solver   = std::shared_ptr<Solver>(),
						 const Parameters params                         = Parameters() ) :
		linear_solver_(linear_solver), active_set_tol_(1e-15)
		{
			set_parameters(params);
		}
		
		bool solve(Vector &x, const Matrix &A, const Vector &b, const Vector &g)
		{
			std::cerr << "[Warning][Deprecated] SemismoothNewton: use the new box constraint interface. This method will be removed shortly" << std::endl;
			std::cout << "[Warning][Deprecated] SemismoothNewton: use the new box constraint interface. This method will be removed shortly" << std::endl;
			
			set_box_constraints(make_upper_bound_constraints(std::make_shared<Vector>(g)));
			return solve(A, b, x);
		}

		inline void set_active_set_tol(const Scalar tol)
		{
			active_set_tol_ = tol;
		}
		
		
		bool solve(const Matrix &A, const Vector &b, Vector &x)  override
		{
			if( constraints_.has_upper_bound() && constraints_.has_lower_bound()) {
				box_solve(A, b, x);
			} else if((constraints_.has_upper_bound() && !constraints_.has_lower_bound()) || (!constraints_.has_upper_bound() && constraints_.has_lower_bound())) {
				single_bound_solve(A, b, x);
			} else {
				std::cout<<"if you do not have constraints, use something else..... \n";
			}
			
			return true;
		}
		
		virtual void set_parameters(const Parameters params) override
		{
			IterativeSolver<Matrix, Vector>::set_parameters(params);
		}
		
		virtual bool set_box_constraints(const BoxConstraints & box)
		{
			constraints_ = box;
			return true;
		}
		
		virtual bool  get_box_constraints(BoxConstraints & box)
		{
			box = constraints_;
			return true;
		}
		
	private:
		// We separate cases with 1 and 2 constraints in order to avoid usless computations in single constraint case
		bool single_bound_solve(const Matrix &A, const Vector &b, Vector &x_new)
		{			
			bool is_upper_bound = constraints_.has_upper_bound();
			
			Vector g;
			if(is_upper_bound) {
				g = *constraints_.upper_bound();
			} else {
				g = *constraints_.lower_bound();
			}
			
			const SizeType local_N = local_size(x_new).get(0);
			
			SizeType iterations = 0;
			bool converged = false;
			
			Vector lambda = local_zeros(local_N);
			Vector active = local_zeros(local_N);
			
			Vector x_old = x_new;
			Vector d, prev_active;
			
			// active/inactive constraints
			Matrix A_c;
			Matrix I_c;
			
			if(is_sparse<Matrix>::value) {
				A_c = local_sparse(local_N, local_N, 1);
				I_c = local_sparse(local_N, local_N, 1);
			} else {
				A_c = local_zeros({local_N, local_N});
				I_c = local_zeros({local_N, local_N});
			}
			
			Scalar f_norm = 9e9;
			
			if(this->verbose())
				this->init_solver("SEMISMOOTH NEWTON METHOD", {" it. ", "|| g ||"});
			
			while(!converged)
			{
				d = lambda + (x_new - g);
				
				{
					Write<Vector> w_a(active);
					Write<Matrix> w_A_c(A_c);
					Write<Matrix> w_I_c(I_c);
					
					Read<Vector> r_d(d);
					
					const Range rr = row_range(A_c);
					
					if(is_upper_bound) {
						for (SizeType i = rr.begin(); i != rr.end(); i++) {
							if (d.get(i) >= -active_set_tol_) {
								A_c.set(i, i, 1.0);
								active.set(i, 1.0);
								
								I_c.set(i, i, 0.0);
							} else {
								I_c.set(i, i, 1.0);
								active.set(i, 0.0);
								
								A_c.set(i, i, 0.0);
							}
						}
					} else {
						//is_lower_bound
						for (SizeType i = rr.begin(); i != rr.end(); i++) {
							if (d.get(i) <= active_set_tol_) {
								A_c.set(i, i, 1.0);
								active.set(i, 1.0);
								
								I_c.set(i, i, 0.0);
							} else {
								I_c.set(i, i, 1.0);
								active.set(i, 0.0);
								
								A_c.set(i, i, 0.0);
							}
						}
					}
				}
				
				if (iterations > 0) {
					const SizeType n_changed = Scalar(norm1(prev_active - active));
					
					if (n_changed == 0) {
						// active set doesn't change anymore => converged
						// fix this to be done in other way
						converged = this->check_convergence(iterations, 1e-15, 1, 1);
						return true;
					}
				}
				
				prev_active = active;
				
				Matrix H = A_c + I_c * A;
				Vector sub_g = (I_c * b + A_c * g);

				assert(!has_nan_or_inf(H));
				assert(!has_nan_or_inf(g));
				
				if(!linear_solver_->solve(H, sub_g, x_new))
					return false;

				if(has_nan_or_inf(x_new)) {
					write("H.m", H);
					write("g.m", sub_g);
					assert(!has_nan_or_inf(x_new));
				}

				lambda = A_c * (b - A * x_new);
				
				assert(!has_nan_or_inf(lambda));

				f_norm = norm2(x_new - x_old);
				
				// print iteration status on every iteration
				if(this->verbose())
					PrintInfo::print_iter_status(iterations, {f_norm});
				
				converged = this->check_convergence(iterations, f_norm, 1, 1);
				
				x_old = x_new;
				iterations++;
			}
			
			return true;
		}
		
		bool box_solve(const Matrix &A, const Vector &b, Vector &x)
		{
			using namespace utopia;
			
			const Size s_A = local_size(A);
			const SizeType n = s_A.get(0);
			const SizeType m = s_A.get(1);
			Scalar x_diff_norm = 0.0;
			
			SizeType it = 0;
			bool converged = false;
			
			const Vector &lb = *constraints_.lower_bound();
			const Vector &ub = *constraints_.upper_bound();
			
			Vector lambda_p = local_zeros(n);
			Vector lambda_m = local_zeros(n);
			
			Vector active_m = local_zeros(n);
			Vector active_p = local_zeros(n);
			
			Vector x_old = x;
			Vector d_p, d_m;
			Vector prev_active_p, prev_active_m;
			Vector active;
			
			// active/inactive constraints
			Matrix A_c_p, A_c_m, A_s, I_c;
			Matrix H;
			
			if(this->verbose()) {
				this->init_solver("SEMISMOOTH NEWTON METHOD", {" it. ", "      || x_k - x_{k-1} ||"});
			}
			
			if(is_sparse<Matrix>::value) {
				A_c_p = 0. * local_identity(n, m);
				A_c_m = 0. * local_identity(n, m);
				I_c = local_identity(n, m);
			} else {
				A_c_p = local_zeros(n);
				A_c_m = local_zeros(n);
				I_c   = local_identity(n, m);
			}
			
			while(!converged) {
				d_p = lambda_p + x - ub;
				d_m = lambda_m + x - lb;
				
				{
					Write<Matrix> w_ac_p(A_c_p);
					Write<Matrix> w_ac_m(A_c_m);
					Write<Vector> w_a_m(active_m);
					Write<Vector> w_a_p(active_p);
					
					Read<Vector> rdp(d_p);
					Read<Vector> rdm(d_m);
					Read<Vector> rb(b);
					
					Range rr = row_range(A_c_p);
					for (SizeType i = rr.begin(); i != rr.end(); i++)
					{
						if (d_p.get(i) >= -active_set_tol_) {
							A_c_p.set(i, i, 1.0);
							
							active_p.set(i, 1.0);
							active_m.set(i, 0.0);
							
						} else if(d_m.get(i) <= active_set_tol_) {
							A_c_m.set(i, i, 1.0);
							
							active_m.set(i, 1.0);
							active_p.set(i, 0.0);
							
						} else {
							A_c_m.set(i, i, 0.0);
							A_c_p.set(i, i, 0.0);
							
							active_m.set(i, 0.0);
							active_p.set(i, 0.0);
						}
					}
				}
				
				active = active_m + active_p;
				A_s = diag(active);
				I_c = local_identity(n, m);
				I_c -= A_s;
				
				if (it > 0) {
					// active set doesn't change anymore => converged
					const SizeType lb_changed = Scalar(norm1(prev_active_m - active_m));
					const SizeType ub_changed = Scalar(norm1(prev_active_p - active_p));
					const SizeType n_changed = lb_changed + ub_changed;
					
					if (n_changed == 0) {
						if(this->verbose() && mpi_world_rank() == 0) {
							std::cout<<"SemismoothNewton:: set is not changing ... \n";
						}
						
						converged = true;
						return true;
					}
				}
				
				prev_active_m = active_m;
				prev_active_p = active_p;
				
				H = A_s + I_c * A;
				
				Vector rhs =  I_c * b + A_c_p * ub + A_c_m * lb;
				linear_solver_->solve(H, rhs, x);
				
				lambda_p = A_c_p * (b - A * x);
				lambda_m = A_c_m * (b - A * x);
				
				// check for change in iterates
				x_diff_norm = norm2(x - x_old);
				
				// print iteration status on every iteration
				if(this->verbose()) {
					PrintInfo::print_iter_status(it, {x_diff_norm});
				}
				
				if(!converged) {
					converged = this->check_convergence(it, x_diff_norm, 1, 1);
				}
				
				x_old = x;
				it++;
			}
			
			return true;
		}
		
		std::shared_ptr <Solver>        linear_solver_;
		BoxConstraints                  constraints_;
		Scalar active_set_tol_;
	};
	
}

#endif //UTOPIA_SOLVER_SEMISMOOTH_NEWTON_HPP
