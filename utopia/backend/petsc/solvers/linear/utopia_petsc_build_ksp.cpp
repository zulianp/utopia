/*
* @Author: kopanicakova
* @Date:   2018-03-08 22:57:34
* @Last Modified by:   kopanicakova
* @Last Modified time: 2018-03-13 22:11:02
*/


#include "utopia_petsc_build_ksp.hpp"
#include "utopia_Instance.hpp"

#undef __FUNCT__
#define __FUNCT__ "KSPSetUp_UTOPIA"
static PetscErrorCode KSPSetUp_UTOPIA(KSP ksp)
{
  PetscFunctionBegin;

  KSP_UTOPIA         *utopia_Pt = (KSP_UTOPIA*)ksp->data;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "KSPSolve_UTOPIA"
static PetscErrorCode KSPSolve_UTOPIA(KSP ksp)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
 
  Mat            Amat, Pmat;

  PCGetOperators(ksp->pc, &Amat, &Pmat);

  KSP_UTOPIA *utopia_ls = (KSP_UTOPIA*)ksp->data;
  
  utopia_ls->utopia_set_tolerances(ksp->rtol, ksp->abstol, ksp->divtol, ksp->max_it); 
  utopia_ls->utopia_solve_routine(Amat, Pmat,  ksp->vec_rhs, ksp->vec_sol); 

  // petsc preconditioner combined with utopia solver ... 
  utopia_ls->get_convergence_reason(ksp->its, ksp->reason); 

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KSPSetSolveRoutine_UTOPIA"
PetscErrorCode  KSPSetSolveRoutine_UTOPIA(KSP ksp, std::function< void(const Mat &, const Mat &, const Vec &, Vec &) > solve_routine)
{
  PetscFunctionBegin;

  KSP_UTOPIA *utopia_ls = (KSP_UTOPIA*)ksp->data;
  utopia_ls->utopia_solve_routine = solve_routine;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KSPSetTolerances_UTOPIA"
PetscErrorCode KSPSetTolerances_UTOPIA(KSP ksp, std::function< void(const PetscReal &, const PetscReal &, const PetscReal &, const PetscInt &) > set_tolerances)
{
  PetscFunctionBegin;

  KSP_UTOPIA *utopia_ls = (KSP_UTOPIA*)ksp->data;
  utopia_ls->utopia_set_tolerances = set_tolerances;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KSPSetGetConvergenceReason_UTOPIA"
PetscErrorCode KSPSetGetConvergenceReason_UTOPIA(KSP ksp, std::function< void(PetscInt &,  KSPConvergedReason&) > convergence_reason)
{
  PetscFunctionBegin;

  KSP_UTOPIA *utopia_ls = (KSP_UTOPIA*)ksp->data;
  utopia_ls->get_convergence_reason = convergence_reason;

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "KSPDestroy_UTOPIA"
PetscErrorCode KSPDestroy_UTOPIA(KSP ksp)
{
  PetscFunctionBegin;
  KSP_UTOPIA         *utopia = (KSP_UTOPIA*)ksp->data;

  KSPDestroyDefault(ksp);

  PetscObjectComposeFunction((PetscObject)ksp, "KSPUTOPIASetSolveRoutine_C",NULL);
  PetscObjectComposeFunction((PetscObject)ksp, "KSPSetTolerancesRoutine_C",NULL);
  PetscObjectComposeFunction((PetscObject)ksp, "KSPSetGetConvergenceReason_C",NULL);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KSPView_UTOPIA"
PetscErrorCode KSPView_UTOPIA(KSP ksp,PetscViewer viewer)
{

  PetscFunctionBegin;


  PetscBool      iascii;
  PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);
  KSP_UTOPIA *utopia_ls = (KSP_UTOPIA*)ksp->data;

  // maybe som future printouts ... 
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KSPSetFromOptions_UTOPIA"
PetscErrorCode KSPSetFromOptions_UTOPIA(PetscOptionItems *PetscOptionsObject,KSP ksp)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;

  // not much to add so far... 
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "KSPCreate_UTOPIA"
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

  ierr = PetscObjectComposeFunction((PetscObject)ksp, "KSPUTOPIASetSolveRoutine_C", KSPSetSolveRoutine_UTOPIA);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)ksp, "KSPSetTolerancesRoutine_C", KSPSetTolerances_UTOPIA);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)ksp, "KSPSetGetConvergenceReason_C", KSPSetGetConvergenceReason_UTOPIA);CHKERRQ(ierr);

  m_utopia_warning_once("> FIXME: setting KSPSetSupportedNorm(ksp, KSP_NORM_NATURAL, PC_RIGHT, 1) see:\n"
                        "  http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetSupportedNorm.html");
    
  ierr = KSPSetSupportedNorm(ksp, KSP_NORM_NATURAL, PC_RIGHT, 1); CHKERRQ(ierr);

  //Maybe use this?
  // ierr = KSPSetSupportedNorm(ksp, KSP_NORM_UNPRECONDITIONED, 1); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}