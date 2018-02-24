/*
 20180222 CP

 @ PZ: stuff you need to set:
 -----------------------------
 1: call NBGSCreate(...)
 2: call NBGSStep(...)
 3: call NBGSDestroy(...)

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

    PetscInt m, M, n, N, size, rank;

    Mat * Dinv;
    Mat Aii; /* system matrix  */
    Vec f, fi; /* the right-hand side */

    Mat loc_A; // holds the local diagonal part of A
    Vec x, res, u, lb_u, ub_u, loc_res, loc_x, loc_u, loc_lb, loc_ub, c, Ac, mask;

    PetscBool constrained;
    PetscBool with_ls;

    PetscInt blocksize;
    PetscScalar alpha, num_active_coordinates, num_changed_coordinates;

    Vec active_coordinates;  // contains the indivdual coordinate indices of dirichlet nodes (also includes the nodes in contact)
    Vec prev_active_coordinates;
    Vec changed_coordinates;

    IS is_nodes;
    IS is_inactive_nodes;

    Vec ub = nullptr;
    Vec lb = nullptr;

    PetscErrorCode (*constrain)(MatScalar * diag, PetscScalar * t, PetscScalar * x,
                                PetscScalar * lb, PetscScalar * ub, int row, int bs, PetscBool * constrained);

    PetscBool initialized;

} NBGS_CTX;

#endif

extern PetscErrorCode NBGSDestroy(NBGS_CTX * nbgs);
extern PetscErrorCode NBGSCreate(NBGS_CTX * nbgs, Mat _A, Vec _x, Vec _ub, Vec _lb, PetscInt _blocksize, PetscBool _with_ls, PetscBool _constrained, PetscErrorCode (*constrain)(MatScalar * diag, PetscScalar * t, PetscScalar * x,PetscScalar * lb, PetscScalar * ub, int row, int bs, PetscBool * constrained));
extern PetscErrorCode NBGSStep(NBGS_CTX * nbgs, Mat A, Vec b, Vec x, Vec lb, Vec ub ,PetscInt _blocksize);



