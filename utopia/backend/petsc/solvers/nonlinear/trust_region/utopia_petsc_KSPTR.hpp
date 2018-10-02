#ifndef UTOPIA_TR_SUBPROBLEM_KSP_TR_HPP
#define UTOPIA_TR_SUBPROBLEM_KSP_TR_HPP
#include "utopia_TRSubproblem.hpp"
#include <string>

namespace utopia
{
    template<typename Matrix, typename Vector>
    class SteihaugToint<Matrix, Vector, PETSC>  : 	virtual public TRSubproblem<Matrix, Vector>,
    			   									virtual public KSPSolver<Matrix, Vector>
	{

    	typedef UTOPIA_SCALAR(Vector) Scalar;
		typedef utopia::KSPSolver<Matrix, Vector> KSPSolver;
		typedef utopia::TRSubproblem<Matrix, Vector> TRSubproblem;
		typedef utopia::Preconditioner<Vector> Preconditioner;
		typedef utopia::IterativeSolver<Matrix, Vector> IterativeSolver;

		static_assert(Traits<Matrix>::Backend == utopia::PETSC, "utopia::KSP_TR:: only works with petsc types");

    public:
        SteihaugToint(const Parameters params = Parameters(), const std::string &preconditioner = "jacobi")
        : TRSubproblem(params), KSPSolver(params)
        {
            this->pc_type(preconditioner);
            this->ksp_type("stcg");
        }

        virtual SteihaugToint<Matrix, Vector, PETSC>* clone() const override {
        	return new SteihaugToint<Matrix, Vector, PETSC>(*this);
        }

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
	  	virtual void set_preconditioner(const std::shared_ptr<Preconditioner> &precond) override
	   	{
	    	KSPSolver::set_preconditioner(precond);
	   	}


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
		    	ierr = KSPSTCGSetRadius(ksp, this->current_radius());

			#else
		    	KSPCGSetRadius(ksp, this->current_radius());
			#endif
	        ierr = KSPSetTolerances(ksp, KSPSolver::rtol(), KSPSolver::atol(), PETSC_DEFAULT,  TRSubproblem::max_it());
	    }

    };

    template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class Nash {};

    template<typename Matrix, typename Vector>
    class Nash<Matrix, Vector, PETSC>  : 	virtual public TRSubproblem<Matrix, Vector>,
    			   								virtual public KSPSolver<Matrix, Vector>
	{

    	typedef UTOPIA_SCALAR(Vector) Scalar;
		typedef utopia::KSPSolver<Matrix, Vector> KSPSolver;
		typedef utopia::TRSubproblem<Matrix, Vector> TRSubproblem;
		typedef utopia::Preconditioner<Vector> Preconditioner;
		typedef utopia::IterativeSolver<Matrix, Vector> IterativeSolver;

		static_assert(Traits<Matrix>::Backend == utopia::PETSC, "utopia::KSP_TR:: only works with petsc types");

    public:
        Nash(const Parameters params = Parameters(), const std::string &preconditioner = "jacobi")
        : TRSubproblem(params), KSPSolver(params)
        {
            this->pc_type(preconditioner);
            this->ksp_type("nash");
        }

        virtual Nash<Matrix, Vector, PETSC>* clone() const override {
        	return new Nash<Matrix, Vector, PETSC>(*this);
        }


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
	  	virtual void set_preconditioner(const std::shared_ptr<Preconditioner> &precond) override
	   	{
	    	KSPSolver::set_preconditioner(precond);
	   	}


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
				ierr = KSPNASHSetRadius(ksp, this->current_radius());

			#else
				KSPCGSetRadius(ksp, this->current_radius());
			#endif
	        ierr = KSPSetTolerances(ksp, KSPSolver::rtol(), KSPSolver::atol(), PETSC_DEFAULT,  TRSubproblem::max_it());
	    }

    };

    template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class Lanczos {};

   	template<typename Matrix, typename Vector>
    class Lanczos<Matrix, Vector, PETSC> : 	virtual public TRSubproblem<Matrix, Vector>,
    			   									virtual public KSPSolver<Matrix, Vector>
	{

    	typedef UTOPIA_SCALAR(Vector) Scalar;
		typedef utopia::KSPSolver<Matrix, Vector> KSPSolver;
		typedef utopia::TRSubproblem<Matrix, Vector> TRSubproblem;
		typedef utopia::Preconditioner<Vector> Preconditioner;
		typedef utopia::IterativeSolver<Matrix, Vector> IterativeSolver;

		static_assert(Traits<Matrix>::Backend == utopia::PETSC, "utopia::KSP_TR:: only works with petsc types");

    public:
        Lanczos(const Parameters params = Parameters(), const std::string &preconditioner = "jacobi")
        : TRSubproblem(params), KSPSolver(params)
        {
            this->pc_type(preconditioner);
            this->ksp_type("nash");
        }


        virtual Lanczos<Matrix, Vector, PETSC>* clone() const override {
        	return new Lanczos<Matrix, Vector, PETSC>(*this);
        }

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
	  	virtual void set_preconditioner(const std::shared_ptr<Preconditioner> &precond) override
	   	{
	    	KSPSolver::set_preconditioner(precond);
	   	}


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
					        ierr = KSPGLTRSetRadius(ksp, this->current_radius());
			#else
					    KSPCGSetRadius(ksp, this->current_radius());
			#endif
	        ierr = KSPSetTolerances(ksp, KSPSolver::rtol(), KSPSolver::atol(), PETSC_DEFAULT,  TRSubproblem::max_it());
	    }

    };


    template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class CGNE {};

    template<typename Matrix, typename Vector>
    class CGNE<Matrix, Vector, PETSC> final : 	virtual public TRSubproblem<Matrix, Vector>,
    			   								virtual public KSPSolver<Matrix, Vector>
	{
		typedef UTOPIA_SCALAR(Vector) Scalar;
		typedef utopia::KSPSolver<Matrix, Vector> KSPSolver;
		typedef utopia::TRSubproblem<Matrix, Vector> TRSubproblem;
		typedef utopia::Preconditioner<Vector> Preconditioner;
		typedef utopia::IterativeSolver<Matrix, Vector> IterativeSolver;

		static_assert(Traits<Matrix>::Backend == utopia::PETSC, "utopia::KSP_TR:: only works with petsc types");


    public:
        CGNE(const Parameters params = Parameters(), const std::string &preconditioner = "jacobi")
        : TRSubproblem(params), KSPSolver(params)
        {
            this->pc_type(preconditioner);
            this->ksp_type("cgne");
        }


        virtual CGNE<Matrix, Vector, PETSC>* clone() const override {
        	return new CGNE<Matrix, Vector, PETSC>(*this);
        }

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
	  	virtual void set_preconditioner(const std::shared_ptr<Preconditioner> &precond) override
	   	{
	    	KSPSolver::set_preconditioner(precond);
	   	}


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
			#else
					    KSPCGSetRadius(ksp, this->current_radius());
			#endif

	        ierr = KSPSetTolerances(ksp, KSPSolver::rtol(), KSPSolver::atol(), PETSC_DEFAULT,  TRSubproblem::max_it());
	    }

    };

}

#endif //UTOPIA_TR_SUBPROBLEM_KSP_TR_HPP
