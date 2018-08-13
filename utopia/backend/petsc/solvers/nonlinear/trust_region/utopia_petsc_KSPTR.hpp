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

    	KSP_TR(const Parameters params = Parameters()): TRSubproblem(params), KSPSolver(params)//,{"stcg", "nash", "cgne", "gltr", "qcg"})
        { 
        	set_parameters(params); 
        	this->ksp_type("stcg");
        	this->pc_type("jacobi");
        }

        KSP_TR(const std::string type, const Parameters params = Parameters()): 
    															TRSubproblem(params),
    															KSPSolver(params)
        { 
        	set_parameters(params); 
        	this->pc_type("jacobi");
        	this->ksp_type(type);
        }


        virtual ~KSP_TR(){}

	    /**
	     * @brief      Sets the parameters.
	     *
	     * @param[in]  params  The parameters
	     */
	    virtual void set_parameters(const Parameters params) override
	    {
	    	Parameters params_copy = params;
	    	params_copy.preconditioner_type(this->pc_type().c_str());
	        KSPSolver::set_parameters(params_copy);
	        TRSubproblem::set_parameters(params_copy);
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


	public:
	    void atol(const Scalar & atol_in)  {  KSPSolver::atol(atol_in); }; 
        void rtol(const Scalar & rtol_in)  {  KSPSolver::rtol(rtol_in); }; 
        void stol(const Scalar & stol_in)  {  KSPSolver::stol(stol_in); }; 	   	
        
        Scalar      atol() const               	{ return KSPSolver::atol(); } 
        Scalar      rtol()  const              	{ return KSPSolver::rtol(); } 
        Scalar      stol()  const              	{ return KSPSolver::stol(); }
        

        virtual KSP_TR<Matrix, Vector, PETSC> * clone() const override {
        	return new KSP_TR(this->ksp_type());
        }
        
	// protected:

        virtual bool apply(const Vector &b, Vector &x) override
    	{
    		Vector grad = -1 * b; 
    		KSPSolver::apply(grad, x); 
    		return true; 
    	}

	    /**
	     * @brief      Update function. 
	     */
	    virtual void update(const std::shared_ptr<const Matrix> &op) override
	    {
	         KSPSolver::update(op);
	         set_ksp_options(KSPSolver::implementation());
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
	    	PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
	    	ierr = KSPSetFromOptions(ksp); 

	        ierr = KSPSetType(ksp, this->ksp_type().c_str()); 
	        ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);

	        if(!this->get_preconditioner()) 
	        {
	            PC pc; 
	            ierr = KSPGetPC(ksp, &pc);
	            ierr = PCSetType(pc, this->pc_type().c_str());
	        }
	        
#if UTOPIA_PETSC_VERSION_LESS_THAN(3,8,0) 
			if(this->ksp_type() == "qcg")
				ierr = KSPQCGSetTrustRegionRadius(ksp, this->current_radius()); 
		    else if(this->ksp_type() == "gltr")    
		        ierr = KSPGLTRSetRadius(ksp, this->current_radius()); 
			else if(this->ksp_type() == "nash")    	        
		        ierr = KSPNASHSetRadius(ksp, this->current_radius()); 
			else
		        ierr = KSPSTCGSetRadius(ksp, this->current_radius()); 

#else
		    KSPCGSetRadius(ksp, this->current_radius());
#endif
	        ierr = KSPSetTolerances(ksp, KSPSolver::rtol(), KSPSolver::atol(), PETSC_DEFAULT,  TRSubproblem::max_it());
	    }

    };

    template<typename Matrix, typename Vector>
    class SteihaugToint<Matrix, Vector, PETSC> : public KSP_TR<Matrix, Vector, PETSC> {
    public:
        SteihaugToint(const Parameters params = Parameters(), const std::string &preconditioner = "jacobi")
        : KSP_TR<Matrix, Vector, PETSC>(params)
        {
            this->pc_type(preconditioner);
            this->ksp_type("stcg");
        }
        
        inline void set_parameters(const Parameters params) override
        {
            Parameters params_copy = params;
            params_copy.lin_solver_type("stcg");
            params_copy.preconditioner_type(this->pc_type().c_str());
            KSP_TR<Matrix, Vector, PETSC>::set_parameters(params_copy);
        }
    };
    
}

#endif //UTOPIA_TR_SUBPROBLEM_KSP_TR_HPP
