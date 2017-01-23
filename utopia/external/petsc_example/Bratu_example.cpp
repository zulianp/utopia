
static char help[] = "Bratu nonlinear PDE in 2d. \n";

/* ------------------------------------------------------------------------

    Solid Fuel Ignition (SFI) problem.  This problem is modeled by
    the partial differential equation

            -Laplacian u - lambda*exp(u) = 0,  0 < x,y < 1,

    with boundary conditions

             u = 0  for  x = 0, x = 1, y = 0, y = 1.

    A finite difference approximation with the usual 5-point stencil
    is used to discretize the boundary value problem to obtain a nonlinear
    system of equations.


    make clean && make -j6 &&  mpirun -n 1 ./nonlinear_ex   -da_grid_x 5 -da_grid_y 5 -par 6.0

  ------------------------------------------------------------------------- */
#include <petscdm.h>
#include <petscdmda.h>
#include <petscsnes.h>
#include <petscmatlab.h>
#include <utopia.hpp>   
#include "PetscUtopiaFunction.hpp"   

/*
   User-defined application context - contains data needed by the
   application-provided call-back routines, FormJacobianLocal() and
   FormFunctionLocal().
*/
typedef struct 
{
  PassiveReal param;          /* test problem parameter */
} AppCtx;

/*
   User-defined routines
*/
extern PetscErrorCode FormInitialGuess(DM,AppCtx*,Vec);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo*,PetscScalar**,PetscScalar**,AppCtx*);
extern PetscErrorCode FormJacobianLocal(DMDALocalInfo*,PetscScalar**,Mat,Mat,AppCtx*);
extern PetscErrorCode FormObjectiveLocal(DMDALocalInfo*,PetscScalar**,PetscReal*,AppCtx*);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  SNES           snes;                         /* nonlinear solver */
  Vec            x;                            /* solution vector */
  AppCtx         user;                         /* user-defined work context */
  PetscInt       its;                          /* iterations for convergence */
  PetscErrorCode ierr;
  PetscReal      bratu_lambda_max = 6.81;
  PetscReal      bratu_lambda_min = 0.;
  PetscBool      flg              = PETSC_FALSE;
  DM             da;
  KSP            ksp;
  PetscInt       lits,slits;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize program
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscInitialize(&argc,&argv,(char*)0,help);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize problem parameters
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  user.param = 6.0;
  ierr       = PetscOptionsGetReal(NULL,"-par",&user.param,NULL);CHKERRQ(ierr);
  if (user.param >= bratu_lambda_max || user.param <= bratu_lambda_min) SETERRQ3(PETSC_COMM_SELF,1,"Lambda, %g, is out of range, [%g, %g]", user.param, bratu_lambda_min, bratu_lambda_max);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create nonlinear solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
  ierr = SNESSetCountersReset(snes,PETSC_FALSE);CHKERRQ(ierr);


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create distributed array (DMDA) to manage parallel grid and vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,-4,-4,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,&da);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);CHKERRQ(ierr);
  ierr = DMSetApplicationContext(da,&user);CHKERRQ(ierr);
  ierr = SNESSetDM(snes,da);CHKERRQ(ierr);
  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Extract global vectors from DMDA; then duplicate for remaining
     vectors that are the same types
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set local function evaluation routine
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMDASNESSetFunctionLocal(da,INSERT_VALUES,(DMDASNESFunction)FormFunctionLocal,&user);CHKERRQ(ierr);
  ierr = DMDASNESSetObjectiveLocal(da,(DMDASNESObjective)FormObjectiveLocal,&user);CHKERRQ(ierr);

  ierr = PetscOptionsGetBool(NULL,"-fd",&flg,NULL);CHKERRQ(ierr);
  if (!flg) 
  {
    ierr = DMDASNESSetJacobianLocal(da,(DMDASNESJacobian)FormJacobianLocal,&user);CHKERRQ(ierr);
  }

  


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Customize nonlinear solver; set runtime options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  ierr = FormInitialGuess(da,&user,x);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Solve nonlinear system
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    printf ("BRATU_EXAMPLE:: solve_test... \n ");

    Vec gradient; 
    Mat jac; 

    VecDuplicate(x, &gradient); 
    SNESComputeFunction(snes, x, gradient);   

    PetscInt dim; 
    VecGetSize(x, &dim); 

    ierr = MatCreate(PETSC_COMM_WORLD,&jac); CHKERRQ(ierr);
    ierr = MatSetSizes(jac, PETSC_DECIDE, PETSC_DECIDE, dim, dim); CHKERRQ(ierr);
    ierr = MatSetFromOptions(jac); CHKERRQ(ierr);
    ierr = MatSetUp(jac); CHKERRQ(ierr);
    ierr = SNESSetJacobian(snes,jac,jac,SNESComputeJacobianDefault,NULL); CHKERRQ(ierr);

    snes->vec_sol_update = x; 
    ierr = SNESComputeJacobian(snes, x, jac, jac);CHKERRQ(ierr);

    {
      using namespace utopia; 
      DSMatrixd jac_utopia; 
      convert(jac, jac_utopia); 
      DVectord x_u; 
      convert(x, x_u); 

      PETSCUtopiaNonlinearFunction<DSMatrixd, DVectord> fun(snes);

      Parameters params; 
      params.tol(1e-8); 
      // params.solver_type("LINE_SEARCH");
      params.solver_type("TRUST_REGION"); 
      params.lin_solver_type(BICGSTAB_TAG); 
      params.trust_region_alg(DOGLEG_TAG);
      params.verbose(true);  
      
      solve(fun, x_u, params); 
    }
  
  

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return 0;
}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormInitialGuess"
/*
   FormInitialGuess - Forms initial approximation.

   Input Parameters:
   user - user-defined application context
   X - vector

   Output Parameter:
   X - vector
 */
PetscErrorCode FormInitialGuess(DM da,AppCtx *user,Vec X)
{
  PetscInt       i,j,Mx,My,xs,ys,xm,ym;
  PetscErrorCode ierr;
  PetscReal      lambda,temp1,temp,hx,hy;
  PetscScalar    **x;

  PetscFunctionBeginUser;
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);

  lambda = user->param;
  hx     = 1.0/(PetscReal)(Mx-1);
  hy     = 1.0/(PetscReal)(My-1);
  temp1  = lambda/(lambda + 1.0);

  /*
     Get a pointer to vector data.
       - For default PETSc vectors, VecGetArray() returns a pointer to
         the data array.  Otherwise, the routine is implementation dependent.
       - You MUST call VecRestoreArray() when you no longer need access to
         the array.
  */
  ierr = DMDAVecGetArray(da,X,&x);CHKERRQ(ierr);

  /*
     Get local grid boundaries (for 2-dimensional DMDA):
       xs, ys   - starting grid indices (no ghost points)
       xm, ym   - widths of local grid (no ghost points)

  */
  ierr = DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);CHKERRQ(ierr);

  /*
     Compute initial guess over the locally owned part of the grid
  */
  for (j=ys; j<ys+ym; j++) {
    temp = (PetscReal)(PetscMin(j,My-j-1))*hy;
    for (i=xs; i<xs+xm; i++) {
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        /* boundary conditions are all zero Dirichlet */
        x[j][i] = 0.0;
      } else {
        x[j][i] = temp1*PetscSqrtReal(PetscMin((PetscReal)(PetscMin(i,Mx-i-1))*hx,temp));
      }
    }
  }

  /*
     Restore vector
  */
  ierr = DMDAVecRestoreArray(da,X,&x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormFunctionLocal"
/*
   FormFunctionLocal - Evaluates nonlinear function, F(x) on local process patch


 */
PetscErrorCode FormFunctionLocal(DMDALocalInfo *info,PetscScalar **x,PetscScalar **f,AppCtx *user)
{
  PetscErrorCode ierr;
  PetscInt       i,j;
  PetscReal      lambda,hx,hy,hxdhy,hydhx,sc;
  PetscScalar    u,ue,uw,un,us,uxx,uyy;

  PetscFunctionBeginUser;
  lambda = user->param;
  hx     = 1.0/(PetscReal)(info->mx-1);
  hy     = 1.0/(PetscReal)(info->my-1);
  sc     = hx*hy*lambda;
  hxdhy  = hx/hy;
  hydhx  = hy/hx;
  /*
     Compute function over the locally owned part of the grid
  */
  for (j=info->ys; j<info->ys+info->ym; j++) {
    for (i=info->xs; i<info->xs+info->xm; i++) {
      if (i == 0 || j == 0 || i == info->mx-1 || j == info->my-1) {
        f[j][i] = 2.0*(hydhx+hxdhy)*x[j][i];
      } else {
        u  = x[j][i];
        uw = x[j][i-1];
        ue = x[j][i+1];
        un = x[j-1][i];
        us = x[j+1][i];

        if (i-1 == 0) uw = 0.;
        if (i+1 == info->mx-1) ue = 0.;
        if (j-1 == 0) un = 0.;
        if (j+1 == info->my-1) us = 0.;

        uxx     = (2.0*u - uw - ue)*hydhx;
        uyy     = (2.0*u - un - us)*hxdhy;
        f[j][i] = uxx + uyy - sc*PetscExpScalar(u);
      }
    }
  }
  ierr = PetscLogFlops(11.0*info->ym*info->xm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormObjectiveLocal"
/*
   FormFunctionLocal - Evaluates nonlinear function, F(x) on local process patch


 */
PetscErrorCode FormObjectiveLocal(DMDALocalInfo *info,PetscScalar **x,PetscReal *obj,AppCtx *user)
{
  PetscErrorCode ierr;
  PetscInt       i,j;
  PetscReal      lambda,hx,hy,hxdhy,hydhx,sc,lobj=0;
  PetscScalar    u,ue,uw,un,us,uxux,uyuy;
  MPI_Comm       comm;

  PetscFunctionBeginUser;
  *obj   = 0;
  ierr = PetscObjectGetComm((PetscObject)info,&comm);CHKERRQ(ierr);
  lambda = user->param;
  hx     = 1.0/(PetscReal)(info->mx-1);
  hy     = 1.0/(PetscReal)(info->my-1);
  sc     = hx*hy*lambda;
  hxdhy  = hx/hy;
  hydhx  = hy/hx;
  /*
     Compute function over the locally owned part of the grid
  */
  for (j=info->ys; j<info->ys+info->ym; j++) {
    for (i=info->xs; i<info->xs+info->xm; i++) {
      if (i == 0 || j == 0 || i == info->mx-1 || j == info->my-1) {
        lobj += PetscRealPart((hydhx + hxdhy)*x[j][i]*x[j][i]);
      } else {
        u  = x[j][i];
        uw = x[j][i-1];
        ue = x[j][i+1];
        un = x[j-1][i];
        us = x[j+1][i];

        if (i-1 == 0) uw = 0.;
        if (i+1 == info->mx-1) ue = 0.;
        if (j-1 == 0) un = 0.;
        if (j+1 == info->my-1) us = 0.;

        /* F[u] = 1/2\int_{\omega}\nabla^2u(x)*u(x)*dx */

        uxux = u*(2.*u - ue - uw)*hydhx;
        uyuy = u*(2.*u - un - us)*hxdhy;

        lobj += PetscRealPart(0.5*(uxux + uyuy) - sc*PetscExpScalar(u));
      }
    }
  }
  ierr = PetscLogFlops(12.0*info->ym*info->xm);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&lobj,obj,1,MPIU_REAL,MPIU_SUM,comm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormJacobianLocal"
/*
   FormJacobianLocal - Evaluates Jacobian matrix on local process patch
*/
PetscErrorCode FormJacobianLocal(DMDALocalInfo *info,PetscScalar **x,Mat jac,Mat jacpre,AppCtx *user)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k;
  MatStencil     col[5],row;
  PetscScalar    lambda,v[5],hx,hy,hxdhy,hydhx,sc;

  PetscFunctionBeginUser;
  lambda = user->param;
  hx     = 1.0/(PetscReal)(info->mx-1);
  hy     = 1.0/(PetscReal)(info->my-1);
  sc     = hx*hy*lambda;
  hxdhy  = hx/hy;
  hydhx  = hy/hx;


  /*
     Compute entries for the locally owned part of the Jacobian.
      - Currently, all PETSc parallel matrix formats are partitioned by
        contiguous chunks of rows across the processors.
      - Each processor needs to insert only elements that it owns
        locally (but any non-local elements will be sent to the
        appropriate processor during matrix assembly).
      - Here, we set all entries for a particular row at once.
      - We can set matrix entries either using either
        MatSetValuesLocal() or MatSetValues(), as discussed above.
  */
  for (j=info->ys; j<info->ys+info->ym; j++) {
    for (i=info->xs; i<info->xs+info->xm; i++) {
      row.j = j; row.i = i;
      /* boundary points */
      if (i == 0 || j == 0 || i == info->mx-1 || j == info->my-1) {
        v[0] =  2.0*(hydhx + hxdhy);
        ierr = MatSetValuesStencil(jacpre,1,&row,1,&row,v,INSERT_VALUES);CHKERRQ(ierr);
      } else {
        k = 0;
        /* interior grid points */
        if (j-1 != 0) {
          v[k]     = -hxdhy;
          col[k].j = j - 1; col[k].i = i;
          k++;
        }
        if (i-1 != 0) {
          v[k]     = -hydhx;
          col[k].j = j;     col[k].i = i-1;
          k++;
        }

        v[k] = 2.0*(hydhx + hxdhy) - sc*PetscExpScalar(x[j][i]); col[k].j = row.j; col[k].i = row.i; k++;

        if (i+1 != info->mx-1) {
          v[k]     = -hydhx;
          col[k].j = j;     col[k].i = i+1;
          k++;
        }
        if (j+1 != info->mx-1) {
          v[k]     = -hxdhy;
          col[k].j = j + 1; col[k].i = i;
          k++;
        }
        ierr = MatSetValuesStencil(jacpre,1,&row,k,col,v,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }

  /*
     Assemble matrix, using the 2-step process:
       MatAssemblyBegin(), MatAssemblyEnd().
  */
  ierr = MatAssemblyBegin(jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  /*
     Tell the matrix we will never add a new nonzero location to the
     matrix. If we do, it will generate an error.
  */
  ierr = MatSetOption(jac,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}