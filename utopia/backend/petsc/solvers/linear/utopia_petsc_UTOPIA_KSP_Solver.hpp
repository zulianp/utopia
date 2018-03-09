/*
* @Author: Alena Kopanicakova
* @Date:   2018-03-3
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2018-03-3
*/
#ifndef UTOPIA_UTOPIA_KSP_HPP
#define UTOPIA_UTOPIA_KSP_HPP


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

#endif