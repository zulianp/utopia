#ifndef UTOPIA_PETSC_BUILD_KSP_HPP
#define UTOPIA_PETSC_BUILD_KSP_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_petsc_Types.hpp"
#include "utopia_petsc.hpp"


#define KSPUTOPIA         "utopia"

/*
        Defines the basic KSP object
*/
#include <petsc/private/kspimpl.h>
#include <iostream>

PETSC_EXTERN PetscErrorCode KSPCreate_UTOPIA(KSP);
PETSC_INTERN PetscErrorCode KSPDestroy_UTOPIA(KSP);
PETSC_INTERN PetscErrorCode KSPReset_UTOPIA(KSP);
PETSC_INTERN PetscErrorCode KSPView_UTOPIA(KSP,PetscViewer);
PETSC_INTERN PetscErrorCode KSPSetFromOptions_UTOPIA(PetscOptionItems *PetscOptionsObject,KSP);

PETSC_INTERN PetscErrorCode KSPSetSolveRoutine_UTOPIA(KSP ksp, std::function< void(const Mat &, const Mat &, const Vec &, Vec &) > solve_routine); 
PETSC_INTERN PetscErrorCode KSPSetTolerances_UTOPIA(KSP ksp, std::function< void(const PetscReal &, const PetscReal &, const PetscReal &, const PetscInt &) > utopia_set_tolerances); 
PETSC_INTERN PetscErrorCode KSPSetGetConvergenceReason_UTOPIA(KSP ksp, std::function< void(PetscInt &,  KSPConvergedReason&) > convergence_reason); 

typedef struct 
{
  std::function< void(const Mat &, const Mat &, const Vec &, Vec &) > utopia_solve_routine; 
  std::function< void(const PetscReal &, const PetscReal &, const PetscReal &, const PetscInt &) > utopia_set_tolerances; 
  std::function< void(PetscInt &,  KSPConvergedReason&) > get_convergence_reason; 

} KSP_UTOPIA;


namespace utopia {
	template<typename Matrix, typename Vector>
	void build_ksp(const std::shared_ptr<utopia::LinearSolver<Matrix, Vector> > & lin_solver, KSP &ksp)
	{
		using namespace utopia; 

		// check if our options overwrite this 
		KSPSetFromOptions(ksp); 
		KSPSetType(ksp, KSPUTOPIA); 

		std::function<void(const Mat &, const Mat &, const Vec &, Vec &)> solve_routine = [&lin_solver](const Mat &A, const Mat & P, const Vec &b, Vec & x)
		{
		    Vector  b_ut, x_ut; 

		    // we need to get some better way how to do this
		    utopia::convert(x, x_ut); 
		    utopia::convert(b, b_ut); 

		    const Matrix A_ut = utopia::sparse_mref(A); 

		    lin_solver->solve(A_ut, b_ut, x_ut); 

			convert(x_ut, x); 
		};


	    std::function< void(const PetscReal &, const PetscReal &, 
	                          const PetscReal &, const PetscInt &) > tol_routine = [&lin_solver](const PetscReal & rtol, 
	                                                                              const PetscReal & abstol, 
	                                                                              const PetscReal & dtol, 
	                                                                              const PetscInt & maxits)
	    {
	        if (dynamic_cast<utopia::IterativeSolver<Matrix, Vector>*>(lin_solver.get()) != nullptr)
	        {
	          auto ls = dynamic_cast<utopia::IterativeSolver<Matrix, Vector> *>(lin_solver.get()); 
	          ls->atol(abstol); 
	          ls->rtol(rtol); 
	          ls->max_it(maxits); 
	        }
	        else
	        {
	          std::cout<<"set tol is not configured for direct solvers yet... \n"; 
	        }
	    };

	    std::function< void(PetscInt &,  KSPConvergedReason&) > convergence_routine = [&lin_solver](PetscInt & max_it,  KSPConvergedReason & reason)
	    {
	      if (dynamic_cast<utopia::IterativeSolver<Matrix, Vector>*>(lin_solver.get()) != nullptr)
	      {
	        auto ls = dynamic_cast<utopia::IterativeSolver<Matrix, Vector> *>(lin_solver.get()); 
	        max_it = ls->get_num_it(); 

	        switch (ls->get_convergence_reason() )
	        {
	            // sucess
	            case utopia::ConvergenceReason::CONVERGED_FNORM_ABS:
	                reason = KSP_CONVERGED_ATOL;
	                break;
	                
	            case utopia::ConvergenceReason::CONVERGED_FNORM_RELATIVE:
	                reason = KSP_CONVERGED_RTOL;
	                break;
	                
	            case utopia::ConvergenceReason::CONVERGED_SNORM_RELATIVE:
	                reason = KSP_CONVERGED_STEP_LENGTH;
	                break;
	                    
	            // fail
	            case utopia::ConvergenceReason::DIVERGED_MAX_IT :
	                reason = KSP_DIVERGED_ITS;
	                break;
	                
	            default :
	                reason = KSP_CONVERGED_RTOL_NORMAL;
	        }
	      }
	      else
	      {
	        std::cout<<"get convergence reason is not configured for direct solvers yet... \n"; 
	      }
	    }; 

	  	KSPSetSolveRoutine_UTOPIA(ksp, solve_routine); 
	  	KSPSetTolerances_UTOPIA(ksp, tol_routine); 
		KSPSetGetConvergenceReason_UTOPIA(ksp, convergence_routine); 

		{
		    // we could hook up PC to precondition utopia solver in future... 
		    PC pc; 
		    KSPGetPC(ksp, &pc); 
		    PCSetType(pc, "none");
		}

		KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
	}
}


#endif