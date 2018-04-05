/*
* @Author: Alena Kopanicakova
* @Date:   2016-09-01
* @Last Modified by:   kopanicakova
* @Last Modified time: 2018-01-26 
*/
#include "utopia_petsc_KSPSolver.hpp"
#include "utopia_petsc_Types.hpp"

#include <petscksp.h>
#include <petscpc.h>


namespace utopia 
{

  /**
   * @brief      Interface between petsc shell preconditioner and utopia solvers/preconditioners.
   */
  PetscErrorCode UtopiaPCApplyShell(PC pc, Vec x, Vec y)
  {
    Preconditioner<DVectord> *shell; 
    PCShellGetContext(pc, (void**)&shell);

    // TODO: ref would be nicer here 
    DVectord x_u, y_u;
    convert(x, x_u); 
    convert(y, y_u); 

    shell->apply(x_u, y_u); 
    convert(y_u, y); 

    return 0;
  }


  PetscErrorCode MyKSPMonitor(KSP ksp,PetscInt it,PetscReal rnorm,void *ctx)
  {

    UTOPIA_TRACE * ut_log;
    ut_log = (UTOPIA_TRACE *) ctx;


    Vec            x;
    PetscReal     conv_rate=0.0; 

    KSPBuildSolution(ksp,NULL,&x);  

    if(it >1)
    {
      Vec nom, denom;

      VecDuplicate(x, &denom); 
      VecDuplicate(x, &nom);
      
      VecCopy(ut_log->x_k_1, nom);
      VecCopy(x, denom);

      VecAYPX(nom, -1, x); 
      VecAYPX(denom, -1, ut_log->x_k_2); 


      PetscReal next, previous; 
      VecNorm(nom,  NORM_2, &next); 
      VecNorm(denom,  NORM_2, &previous); 

      conv_rate = next/previous; 

      VecDestroy(&nom); 
      VecDestroy(&denom); 
    }

    if(it> 0)
      VecCopy(ut_log->x_k_1, ut_log->x_k_2);

      // if(it > 0)
      VecCopy(x, ut_log->x_k_1);

    if(ut_log->compute_cond_number)
    {
      if(it == 0)
        PetscPrintf(PETSC_COMM_WORLD,"it           ||r||                   rho                  cond. number \n");

      PetscReal emax, emin; 
      KSPComputeExtremeSingularValues(ksp, &emax, &emin); 
      PetscPrintf(PETSC_COMM_WORLD,"%D     %14.12e         %14.12e        %14.12e \n", it, rnorm, conv_rate, std::abs(emax)/std::abs(emin));
    }
    else
    {
      if(it == 0)
        PetscPrintf(PETSC_COMM_WORLD,"it           ||r||                   rho      \n");
      
      PetscPrintf(PETSC_COMM_WORLD,"%D     %14.12e         %14.12e \n", it, rnorm, conv_rate);
    }

    return 0;
  }


}

