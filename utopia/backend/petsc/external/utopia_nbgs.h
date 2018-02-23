/*

 20180222 CP

 @ PZ: stuff you need to set:
 -----------------------------

 nbgs.blocksize = pctx->passo_params.blocksize;
 nbgs.constrained = PETSC_TRUE;
 nbgs.constrain = constrain_blocknd;
 nbgs.with_linesearch = pctx->nbgs_ls_on;
 nbgs.ub = snes->xu;
 nbgs.lb = snes->xl;

 initialize the vector that hold info on the active nodes etc.
 -----------

 ierr = VecDuplicate([vector like x], &nbgs.active_coordinates); CHKERRQ(ierr);
 ierr = VecDuplicate([vector like x], &nbgs.prev_active_coordinates); CHKERRQ(ierr);
 ierr = VecDuplicate([vector like x], &nbgs.changed_coordinates); CHKERRQ(ierr);

 ierr = VecScale(nbgs.active_coordinates, 0); CHKERRQ(ierr);
 ierr = VecScale(nbgs.prev_active_coordinates, 0); CHKERRQ(ierr);
 ierr = VecScale(nbgs.changed_coordinates, 0); CHKERRQ(ierr);



  - The routine assumes that your sparse petsc matrix has itself the same blocksize as the one you set.

 have fun, and ask me if there are weird errors. There will be another major version change once I add friction.

 */


#ifndef nbgs_h
#define nbgs_h

#include <petscvec.h>
#include <petscmat.h>
#include <petscsnes.h>

#define  NBGS_TOL 1e-18

typedef struct
{
    
    Mat * Dinv;
    Mat A, Aii; /* system matrix  */
    Vec f, fi; /* the right-hand side */
    PetscBool initialized;
    PetscBool initialized_ls;
    PetscBool constrained;
    PetscBool with_linesearch;
    
    PetscInt blocksize; /*!< in case of using NBGS as smoother */
    PetscScalar alpha;  /*!< dampening or relaxation parameter */
    
    Vec active_coordinates;  // contains the indivdual coordinate indices of dirichlet nodes (also includes the nodes in contact)
    Vec prev_active_coordinates;
    Vec changed_coordinates;

    IS is_nodes;
    IS is_inactive_nodes;
    
    Vec ub = nullptr;
    Vec lb = nullptr;
    
    PetscErrorCode (*constrain)(MatScalar * diag, PetscScalar * t, PetscScalar * x,
                                PetscScalar * lb, PetscScalar * ub, int row, int bs, PetscBool * constrained);
    
    PetscBool is_setup = PETSC_FALSE;
    
} NBGS_CTX;

extern PetscErrorCode NBGS(NBGS_CTX * nbgs, Mat A, Vec b, Vec x, Vec lb, Vec ub ,PetscInt _blocksize);

#endif
