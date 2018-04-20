/*
* @Author: alenakopanicakova
* @Date:   2016-06-06
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-11-06
*/
#ifndef UTOPIA_CONJUGATE_GRAD_H
#define UTOPIA_CONJUGATE_GRAD_H

#include "utopia_IterativeSolver.hpp"
#include "utopia_Parameters.hpp"
#include "utopia_Preconditioner.hpp"
#include "utopia_Smoother.hpp"

#include <memory>

namespace utopia 
{

    /**
     * @brief      Conjugate Gradient solver. Works with all utopia tensor types. 
     * @tparam     Matrix  
     * @tparam     Vector  
     */
	template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
	class ConjugateGradient : public IterativeSolver<Matrix, Vector>, public Smoother<Matrix, Vector>
	{
		typedef UTOPIA_SCALAR(Vector) 	 Scalar;
		typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::LinearSolver<Matrix, Vector> Solver;
        typedef utopia::Preconditioner<Vector> Preconditioner;

	public:

		ConjugateGradient(const Parameters params = Parameters())
        {
            set_parameters(params); 
        }

        /**
         * @brief      Sets the parameters.
         *
         * @param[in]  params  The parameters
         */
        void set_parameters(const Parameters params) override
        {
            IterativeSolver<Matrix, Vector>::set_parameters(params); 
        }

        /**
         * @brief      Solution routine for CG. 
         *
         * @param[in]  b     The right hand side. 
         * @param      x     The initial guess/solution.
         *
         * @return true if the linear system has been solved up to required tollerance. False otherwise
         */
         bool apply(const Vector &b, Vector &x) override
         {
            if(precond_) {
                return preconditioned_solve(*this->get_operator(), b, x);
            } else {
                return unpreconditioned_solve(*this->get_operator(), b, x);
            }
         }

         /**
          * @brief      Sets the preconditioner.
          *
          * @param[in]  precond  The preconditioner.
          */
         void set_preconditioner(const std::shared_ptr<Preconditioner> &precond)
         {
             precond_ = precond;
         }

         /*! @brief if overriden the subclass has to also call this one first
          */
         virtual void update(const std::shared_ptr<const Matrix> &op) override
         {
             IterativeSolver<Matrix, Vector>::update(op);

             if(precond_) {
                auto ls_ptr = dynamic_cast<LinearSolver<Matrix, Vector> *>(precond_.get());
                if(ls_ptr) {
                    ls_ptr->update(op);    
                }
             }
         }

         bool smooth(const Vector &rhs, Vector &x) override
         {
            SizeType temp = this->max_it();
            this->max_it(this->sweeps());
            unpreconditioned_solve(*this->get_operator(), rhs, x);
            this->max_it(temp);
            return true;
         }

         ConjugateGradient * clone() const override
         {
            return new ConjugateGradient(*this);
         }

     private:
        bool unpreconditioned_solve(const Matrix &A, const Vector &b, Vector &x)
        {
        	Scalar it = 0; 
        	Scalar rho = 1., rho_1 = 1., beta = 0., alpha = 1., r_norm = 9e9; 
        	Vector r, p, q; 

            assert(!empty(b));

        	if(empty(x) || size(x).get(0) != size(b).get(0)) {
                x = local_zeros(local_size(b));
                r = b;
            } else {
                assert(local_size(x).get(0) == local_size(b).get(0));
                r = b - A * x; 
            }

            this->init_solver("Utopia Conjugate Gradient", {"it. ", "||r||" }); 
            bool converged = false; 

        	while(!converged)
        	{
        		rho = dot(r, r); 
        		if(it > 0)
        		{
        			beta = rho/rho_1;
        			p = r + beta * p; 
        		}
        		else
        		{   
        			p = r; 
        		}

        		q = A * p; 
        		alpha = rho / dot(p, q); 

        		x += alpha * p; 
        		r -= alpha * q; 

        		rho_1 = rho; 
        		it++; 
        		r_norm = norm2(r); 

                if(this->verbose())
                    PrintInfo::print_iter_status({it, r_norm}); 

                converged = this->check_convergence(it, r_norm, 1, 1); 
                it++; 
        	}

        	return converged; 
        }

        bool preconditioned_solve(const Matrix &A, const Vector &b, Vector &x)
        {
            Scalar it = 0; 
            Scalar beta = 0., alpha = 1., r_norm = 9e9; 
            Vector r, p, Ap, r_new; 

            Vector z     = local_zeros(local_size(b)); 
            Vector z_new = local_zeros(local_size(b)); 

            if(empty(x) || size(x).get(0) != size(b).get(0)) {
                x = local_zeros(local_size(b));
                r = b;
            } else {
                assert(local_size(x).get(0) == local_size(b).get(0));
                r = b - A * x; 
            }

            precond_->apply(r, z);
            p = z; 

            this->init_solver("Utopia Conjugate Gradient", {"it. ", "||r||" }); 
            bool stop = false; 

            while(!stop)
            {
                Ap = A*p;
                alpha = dot(r, z)/dot(p, Ap);
                x += alpha * p;
                r_new = r - alpha * Ap;

                r_norm = norm2(r_new);
                if(r_norm < this->atol()) {
                    if(this->verbose())
                        PrintInfo::print_iter_status({it, r_norm}); 
                    
                    stop = this->check_convergence(it, r_norm, 1, 1); 
                    break;
                }

                precond_->apply(r_new, z_new);
                beta = dot(z_new, r_new)/dot(z, r);

                p = z_new + beta * p;
                r = r_new;
                z = z_new;

                if(this->verbose())
                    PrintInfo::print_iter_status({it, r_norm}); 

                stop = this->check_convergence(it, r_norm, 1, 1); 
                it++; 
            }

            return r_norm <= this->atol(); 
        }


            std::shared_ptr<Preconditioner> precond_;

    };
}


#endif //UTOPIA_CONJUGATE_GRAD_H

