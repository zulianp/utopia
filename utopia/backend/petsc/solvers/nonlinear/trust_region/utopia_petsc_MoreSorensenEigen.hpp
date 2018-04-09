#ifndef UTOPIA_TR_MORE_SORENSEN_EIGEN_HPP
#define UTOPIA_TR_MORE_SORENSEN_EIGEN_HPP
#include "utopia_TRSubproblem.hpp"



namespace utopia 
{

	template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend> 
    class MoreSorensenEigen {};


	/**
	 * @brief      Class for More Sorensen minimization algorithm. Initialization of lambda_0 is based on eigensolution 
	 */
	template<class Matrix, class Vector>
    class MoreSorensenEigen<Matrix, Vector, PETSC> : public TRSubproblem<Matrix, Vector>
    {
		typedef UTOPIA_SCALAR(Vector) Scalar;
		typedef utopia::LinearSolver<Matrix, Vector> 		LinearSolver;
		typedef utopia::EigenValueSlover<Matrix, Vector, PETSC_EXPERIMENTAL> 	EigenSolver;


    public:

    	MoreSorensenEigen(	const std::shared_ptr<LinearSolver> &linear_solver = std::shared_ptr<LinearSolver>(), 
    						const std::shared_ptr<EigenSolver> & eigen_solver = std::shared_ptr<EigenSolver>(), 
    						const Parameters params = Parameters()): 
    						TRSubproblem<Matrix, Vector>(params), 
    						linear_solver_(linear_solver), 
    						eigen_solver_(eigen_solver), 
    						kappa_easy_(0.0001), 
    						max_it_(1000), 
    						lambda_eps_(1e-5)
        {  };

        virtual ~MoreSorensenEigen(){}

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
	        eigen_solver_->solver_type("arnoldi");
	        eigen_solver_->tol(1e-14); 
	        eigen_solver_->solve(H); 

	        eigen_solver_->get_real_eigenpair(0, lambda, eigenvector); 


        	// 	decide if PD case or not
	        lambda = (lambda > 0.0) ? 0.0 : - 1.0 * (lambda - lambda_eps_); 


	        linear_solver_->solve(H, -1 * g, s_k); 
	        s_norm = norm2(s_k); 


	        if(s_norm <= this->current_radius())
	        {

	        	if(lambda == 0.0 || s_norm == this->current_radius())
	        		return true; 
	        	else
	        	{
	        		std::cerr<<"we can not handle hard case  yet ......... \n"; 	
	        		return false; 
	        	}
	        }


	        Matrix H_lambda = H; 

	        for(auto it = 0; it < max_it_; it++)
	        {

		        if( std::abs(s_norm - this->current_radius()) <= kappa_easy_ *  this->current_radius())
		        {
		        	// std::cout<<"we found approximate minimizer after: "<< it << "  iterations. \n"; 
		        	return true; 
		        }

	        	Vector grad_s_lambda = 0 * s_k; 
	        	linear_solver_->solve(H_lambda, -1 * s_k, grad_s_lambda); 
	        	
	        	Scalar grad 	= 1.0/s_norm - 1.0/ this->current_radius(); 
	        	Scalar hessian 	= -1.0 * dot(s_k, grad_s_lambda)/ std::pow(s_norm, 3); 

	        	// new iterate 
	        	lambda -=  grad/hessian; 
		        it++; 

		        // H should be H + \lambda I 
		        // TODO:: investigate why it does not work with mat multiply... 
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