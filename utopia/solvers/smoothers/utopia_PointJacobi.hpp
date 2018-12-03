#ifndef UTOPIA_JACOBI_HPP
#define UTOPIA_JACOBI_HPP
#include "utopia_IterativeSolver.hpp"
#include "utopia_Smoother.hpp"

namespace utopia {

    /**
     * @brief Good example on how to implement a Linear Solver, Preconditioner, and Smoother 
     * within utopia. Precomputations are done in the update method, and, in order to avoid 
     * costly allocations, possible temporaries are stored as member variables.
     */
	template<class Matrix, class Vector>
	class PointJacobi final: public Smoother<Matrix, Vector>, public IterativeSolver<Matrix, Vector> {
		typedef UTOPIA_SCALAR(Vector)    Scalar;
		typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
		typedef utopia::IterativeSolver<Matrix, Vector> Solver;
		typedef utopia::Smoother<Matrix, Vector> Smoother;
		
	public:
		/**
		 * @brief      Very simple Jacobi solver.
		 *
		 * @param[in]  omega  The relaxation parameter (unused atm).
		 */
		PointJacobi(const Parameters params = Parameters())
		{
			set_parameters(params);
		}
		
		bool apply(const Vector &rhs, Vector &x) override
		{
			const Matrix &A = *this->get_operator();
			
			SizeType it = 0;
            r_ = rhs - A * x;
			Scalar g_norm0 = norm2(r_);
			Scalar g_norm = g_norm0;
			SizeType compute_norm_each = 50;
			
			this->init_solver("Point Jacobi", {"it. ", "||r||" });
			
			while(true) {
				sweep(rhs, x);
				
				if(it++ % compute_norm_each == 0) {
                    r_ = rhs - A * x;
					g_norm = norm2(r_);
					
					if(this->verbose()) {
						PrintInfo::print_iter_status(it, {g_norm});
					}
					
					if(this->check_convergence(it, g_norm, g_norm/g_norm0, 1)) {
						return true;
					}
				}
			}
			
			return false;
		}
		
		bool smooth(const Vector &rhs, Vector &x) override
		{			
			for(SizeType it = 0; it < this->sweeps(); it++)
			{
				sweep(rhs, x);
			}
			
			return true;
		}
		
		void set_parameters(const Parameters params) override
		{
			Smoother::set_parameters(params);
			Solver::set_parameters(params);
		}
		
		inline PointJacobi * clone() const override
		{
			return new PointJacobi(*this);
		}
		
		void update(const std::shared_ptr<const Matrix> &op) override
		{
			Solver::update(op);
			
			const auto &A = *op;
			Vector diag_A = diag(A);
			d_inv_ = 1. / diag_A;
			
			// prevents system from being indefinite
			check_indef(d_inv_);
			
			// lower and upper part of A
			LU_ = A;
			LU_ -= Matrix(diag(diag_A));
		}
		
	private:
		Vector d_inv_;
		Matrix LU_;
        Vector r_;
		
		/**
		 * @brief      Checks if there is a zero in the vector, if yes turn it into 1.
		 *
		 * @param      diag_A  { D_{-1}}
		 *
		 * @return     {  }
		 */
		inline bool check_indef(Vector &diag_A)
		{
			ReadAndWrite<Vector> w(diag_A);
			Range rr = range(diag_A);
			
			for (SizeType i = rr.begin(); i != rr.end(); i++)
			{
				if(diag_A.get(i) == 0)
				{
					diag_A.set(i, 1);
				}
			}

			return true;
		}
		
	
		inline bool sweep(const Vector &rhs, Vector &x)
		{
			r_ = rhs - (LU_ * x);
			x = e_mul(d_inv_, r_);
			return true;
		}
		
	};
	
}

#endif //UTOPIA_JACOBI_HPP

