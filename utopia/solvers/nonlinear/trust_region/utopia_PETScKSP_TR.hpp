/*
* @Author: alenakopanicakova
* @Date:   2016-10-04
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-10-13
*/
#ifdef WITH_PETSC

#ifndef UTOPIA_TR_SUBPROBLEM_KSP_TR_HPP
#define UTOPIA_TR_SUBPROBLEM_KSP_TR_HPP
#include "utopia_TRSubproblem.hpp"
#include <string>

namespace utopia 
{

	/**
	 * @brief      Interface to use petsc KSP in TR context \n
	 *				There are 5 different possibilities, grouped as follows: \n
	 *				1. gltr, stcg, nash, qcg: \n
	 										Uses preconditioned conjugate gradient to compute an approximate minimizer of the quadratic function.  \n
											\f$ q(s) = g^T * s + .5 * s^T * H * s   \f$  \n
											subject to  \f$ || s || <= delta  \f$ \n
					2. cgne				: \n 
											Applies the preconditioned conjugate gradient method to the normal equations without explicitly forming A^t*A. \n
										  	Therefore, it's very suitable to use in combination with TR_NormalEquation, LeastSquaresFunction and solution method: dogleg \n
										  	Note, that CGNE is a general-purpose non-symmetric method. \n
	 */
    template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend> 
    class KSP_TR {};


    template<typename Matrix, typename Vector>
    class KSP_TR<Matrix, Vector, PETSC> : 	virtual public TRSubproblem<Matrix, Vector>, 
    			   							virtual public KSPSolver<Matrix, Vector>
    {
		typedef UTOPIA_SCALAR(Vector) Scalar;
		typedef utopia::KSPSolver<Matrix, Vector> KSPSolver;
		typedef utopia::TRSubproblem<Matrix, Vector> TRSubproblem;
		typedef utopia::Preconditioner<Vector> Preconditioner;
		typedef utopia::IterativeSolver<Matrix, Vector> IterativeSolver;

		static_assert(Traits<Matrix>::Backend == utopia::PETSC, "utopia::KSP_TR:: only works with petsc types");

    public:


    	KSP_TR(const Parameters params = Parameters()): 
    															TRSubproblem(params),
    															KSPSolver(params,{"stcg", "nash", "cgne", "gltr", "qcg"})
        { 
        	set_parameters(params); 
        };

        KSP_TR(const std::string type, const Parameters params = Parameters()): 
    															TRSubproblem(params),
    															KSPSolver(params, {type})
        { 
        	set_parameters(params); 
        };


        virtual ~KSP_TR(){}

	    /**
	     * @brief      Sets the parameters.
	     *
	     * @param[in]  params  The parameters
	     */
	    virtual void set_parameters(const Parameters params) override
	    {
	        KSPSolver::set_parameters(params);
	        TRSubproblem::set_parameters(params);
	    }

	   /**
	   * @brief      Sets the preconditioner.
	   *
	   * @param[in]  precond  The precondition
	   */
	  	void set_preconditioner(const std::shared_ptr<Preconditioner> &precond) override
	   	{
	    	KSPSolver::set_preconditioner(precond); 
	   	}


	protected:

        virtual bool apply(const Vector &b, Vector &x) override
    	{
    		KSPSolver::apply(b, x); 
    		x *= -1;  
    		return true; 
    	}

	    /**
	     * @brief      Update function. 
	     */
	    virtual void update(const std::shared_ptr<const Matrix> &op) override
	    {
	         KSPSolver::update(op);
	    }


	private:

	   	/**
	     * @brief      Sets the default options for PETSC KSP solver. \n
	     *             Default: ST-CG
	     *
	     * @param      ksp   The ksp
	     */
	    virtual void set_ksp_options(KSP & ksp) override 
	    {
	        PetscErrorCode ierr;
	        ierr = KSPSetType(ksp, this->ksp_type().c_str()); 
	        ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);

			if(this->ksp_type() == "qcg")
				ierr = KSPQCGSetTrustRegionRadius(ksp, this->current_radius()); 
		    else if(this->ksp_type() == "gltr")    
		        ierr = KSPGLTRSetRadius(ksp, this->current_radius()); 
			else if(this->ksp_type() == "nash")    	        
		        ierr = KSPNASHSetRadius(ksp, this->current_radius()); 
			else
		        ierr = KSPSTCGSetRadius(ksp, this->current_radius()); 

	        ierr = KSPSetTolerances(ksp, KSPSolver::rtol(), KSPSolver::atol(), PETSC_DEFAULT,  KSPSolver::max_it());
	    }

    };

}

#endif //UTOPIA_TR_SUBPROBLEM_KSP_TR_HPP
#endif //WITH_PETSC
