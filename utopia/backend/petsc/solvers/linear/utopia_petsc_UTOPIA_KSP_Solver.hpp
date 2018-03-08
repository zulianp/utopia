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

PETSC_INTERN PetscErrorCode KSPSetSolveRoutine_UTOPIA(KSP ksp, std::function< void(const Mat &, const Vec &, Vec &) > solve_routine); 


typedef struct 
{

  std::function< void(const Mat &, const Vec &, Vec &) > utopia_solve_routine; 

} KSP_UTOPIA;

#endif