#ifndef UTOPIA_TR_MORE_SORENSEN_EIGEN_HPP
#define UTOPIA_TR_MORE_SORENSEN_EIGEN_HPP
#include "utopia_TRSubproblem.hpp"


namespace utopia 
{

	template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend> 
    class MoreSorensenEigen {};

	/**
	 * @brief      Class for More Sorensen minimization algorithm, where initialization of lambda_0 is based on eigen sol. 
	 * 
	 * 			   WARNING:: 	Computation of direction could be more efficient by using cholesky decomposition and store just factors, 
	 * 			   				but it is a bit anoying to do so with petsc => we solve system 2 \times, which is less efficient, 
	 * 			   				but ok as this is just proof of concept solver
	 */
	template<class Matrix, class Vector>
    class MoreSorensenEigen<Matrix, Vector, PETSC> : public TRSubproblem<Matrix, Vector>
    {
		typedef UTOPIA_SCALAR(Vector) Scalar;
		typedef utopia::LinearSolver<Matrix, Vector> 		LinearSolver;
		typedef utopia::EigenSolver<Matrix, Vector, PETSC_EXPERIMENTAL> 	EigenSolver;


    public:
    	MoreSorensenEigen(	const std::shared_ptr<LinearSolver> &linear_solver,  const std::shared_ptr<EigenSolver> & eigen_solver, 
    						const Parameters params = Parameters()): 
    						TRSubproblem<Matrix, Vector>(params), 
    						linear_solver_(linear_solver), 
    						eigen_solver_(eigen_solver), 
    						kappa_easy_(1e-10), 
    						max_it_(1000), 
    						lambda_eps_(1e-5)
        {  };

        virtual ~MoreSorensenEigen(){}

        void kappa_easy(const Scalar & kappa)
        {
        	kappa_easy_ = kappa; 
        }

        Scalar kappa_easy()
        {
        	return kappa_easy_; 
        }

        void lambda_eps(const Scalar & lambda_eps)
        {
        	lambda_eps_ = lambda_eps; 
        }

        Scalar lambda_eps()
        {
        	return lambda_eps_; 
        }

        virtual void set_linear_solver(const std::shared_ptr<LinearSolver > &ls) override
        {
            linear_solver_ = ls; 
        }   

        MoreSorensenEigen * clone() const override
        {
            return new MoreSorensenEigen(std::shared_ptr<LinearSolver>(linear_solver_->clone()), std::shared_ptr<EigenSolver>(eigen_solver_->clone()));
        }


	protected:
        bool unpreconditioned_solve(const Matrix &H, const Vector &g, Vector &s_k) override
        {
        	Scalar lambda, s_norm; 
        	Vector eigenvector; 
        	// init vector... 
        	s_k = 0.0 * g; 

        	// ---------------------- initialization  of lambda_0 ------------------------
        	eigen_solver_->portion_of_spectrum("smallest_real"); 
	        eigen_solver_->number_of_eigenvalues(1); 
	        eigen_solver_->solve(H); 

	        eigen_solver_->get_real_eigenpair(0, lambda, eigenvector); 

        	// 	decide if PD case or not
	        lambda = (lambda > 0.0) ? 0.0 : - 1.0 * (lambda - lambda_eps_); 


	        Matrix H_lambda = H; 
	        if(lambda != 0.0)
	        {
				Write<Matrix> w(H_lambda);
				Range r = row_range(H_lambda);

				for(SizeType i = r.begin(); i != r.end(); ++i) 
					H_lambda.add(i, i, lambda);
			}


	        linear_solver_->solve(H_lambda, -1 * g, s_k); 
	        s_norm = norm2(s_k); 

	        if(s_norm <= this->current_radius())
	        {

	        	if(lambda == 0.0 || s_norm == this->current_radius())
	        		return true;
	        	else
	        	{
	        		// we are in hard case, let's find solution on boundary, which is orthogonal to E_1
	        		//                     because eigenvector is normalized
	        		Scalar alpha = this->quadratic_function(1.0, 2.0 * dot(s_k, eigenvector), dot(s_k, s_k) - (this->current_radius() * this->current_radius())); 
	        		s_k += alpha * eigenvector;  
	        		return true; 
	        	}
	        }

	        for(auto it = 0; it < max_it_; it++)
	        {
		        if( std::abs(s_norm - this->current_radius()) <= kappa_easy_ *  this->current_radius())
		        	return true; 

	        	Vector grad_s_lambda = 0 * s_k; 
	        	linear_solver_->solve(H_lambda, -1 * s_k, grad_s_lambda); 
	        	
	        	Scalar grad 	= 1.0/s_norm - 1.0/ this->current_radius(); 
	        	Scalar hessian 	= -1.0 * dot(s_k, grad_s_lambda)/ std::pow(s_norm, 3); 

	        	
	        	lambda -=  grad/hessian; 

		        // H should be H + \lambda I 
		        // TODO:: investigate why it does not work with mat multiply... something is wrong with mat allocations... 
		        H_lambda = H;
		        {
					Write<Matrix> w(H_lambda);
					Range r = row_range(H_lambda);

			        //You can use set instead of add. [Warning] Petsc does not allow to mix add and set.
					for(SizeType i = r.begin(); i != r.end(); ++i) 
						H_lambda.add(i, i, lambda);
				}

				s_k *= 0.0; 
		       	linear_solver_->solve(H_lambda, -1 * g, s_k); 
	        	s_norm = norm2(s_k); 
		    }

        	return true; 
        }


        bool preconditioned_solve(const Matrix & /*H*/, const Vector &/*g*/, Vector &/*p_k*/) override
        {
        	std::cout<<"MoreSorensenEigen:: preconditioned solve not imlemented yet ... \n"; 
        	return false; 
        }


    private: 
     	std::shared_ptr<LinearSolver> linear_solver_;     
     	std::shared_ptr<EigenSolver> eigen_solver_;     

     	Scalar kappa_easy_; 
     	SizeType max_it_; 
     	Scalar lambda_eps_; 

    };


}

#endif //UTOPIA_TR_MORE_SORENSEN_EIGEN_HPP