/*
* @Author: kopanicakova
* @Date:   2018-03-08 22:57:34
* @Last Modified by:   kopanicakova
* @Last Modified time: 2018-03-09 00:46:04
*/


#include "utopia_petsc_UTOPIA_KSP_Solver.hpp"

static PetscErrorCode KSPSetUp_UTOPIA(KSP ksp)
{
  PetscFunctionBegin;

  KSP_UTOPIA         *utopia_Pt = (KSP_UTOPIA*)ksp->data;
  PetscFunctionReturn(0);
}



static PetscErrorCode KSPSolve_UTOPIA(KSP ksp)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
 
  
  Mat            Amat, Pmat;

  // TODO:: set PC matrix ... 
  PCGetOperators(ksp->pc, &Amat, &Pmat);

  Vec X             = ksp->vec_sol;
  Vec B             = ksp->vec_rhs;

  KSP_UTOPIA *utopia_ls = (KSP_UTOPIA*)ksp->data;
  


  utopia_ls->utopia_solve_routine(Amat, B, X); 



  // TODO:: set-up proper reason ... 
  ksp->reason = KSP_CONVERGED_ATOL;

  PetscFunctionReturn(0);
}


PetscErrorCode  KSPSetSolveRoutine_UTOPIA(KSP ksp, std::function< void(const Mat &, const Vec &, Vec &) > solve_routine)
{
  PetscFunctionBegin;

  KSP_UTOPIA *utopia_ls = (KSP_UTOPIA*)ksp->data;
  utopia_ls->utopia_solve_routine = solve_routine;

  PetscFunctionReturn(0);
}



PetscErrorCode KSPDestroy_UTOPIA(KSP ksp)
{
  PetscFunctionBegin;
  KSP_UTOPIA         *utopia = (KSP_UTOPIA*)ksp->data;

  KSPDestroyDefault(ksp);
  PetscObjectComposeFunction((PetscObject)ksp, "KSPUTOPIASetSolveRoutine_C",NULL);

  PetscFunctionReturn(0);
}


PetscErrorCode KSPView_UTOPIA(KSP ksp,PetscViewer viewer)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}


PetscErrorCode KSPSetFromOptions_UTOPIA(PetscOptionItems *PetscOptionsObject,KSP ksp)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;

  PetscFunctionReturn(0);
}



PETSC_EXTERN PetscErrorCode KSPCreate_UTOPIA(KSP ksp)
{
  PetscErrorCode ierr;
  KSP_UTOPIA         *utopia_solver;

  PetscFunctionBegin;
  ierr = PetscNewLog(ksp, &utopia_solver); CHKERRQ(ierr);

  ksp->data = (void*)utopia_solver;

  ksp->ops->setup          = KSPSetUp_UTOPIA;
  ksp->ops->solve          = KSPSolve_UTOPIA;
  ksp->ops->destroy        = KSPDestroy_UTOPIA;
  ksp->ops->view           = KSPView_UTOPIA;
  ksp->ops->setfromoptions = KSPSetFromOptions_UTOPIA;

  ksp->ops->buildsolution  = KSPBuildSolutionDefault;
  ksp->ops->buildresidual  = KSPBuildResidualDefault;

  ierr = PetscObjectComposeFunction((PetscObject)ksp,"KSPUTOPIASetSolveRoutine_C",KSPSetSolveRoutine_UTOPIA);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}