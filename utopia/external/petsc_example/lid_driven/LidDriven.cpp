
static char help[] = "Nonlinear driven cavity with multigrid in 2d.\n \
  \n\
The 2D driven cavity problem is solved in a velocity-vorticity formulation.\n\
The flow can be driven with the lid or with bouyancy or both:\n\
  -lidvelocity &ltlid&gt, where &ltlid&gt = dimensionless velocity of lid\n\
  -grashof &ltgr&gt, where &ltgr&gt = dimensionless temperature gradent\n\
  -prandtl &ltpr&gt, where &ltpr&gt = dimensionless thermal/momentum diffusity ratio\n\
 -contours : draw contour plots of solution\n\n";
/* in HTML, '&lt' = '<' and '&gt' = '>' */

/*
      See src/ksp/ksp/examples/tutorials/ex45.c
*/

/*T
   Concepts: SNES^solving a system of nonlinear equations (parallel multicomponent example);
   Concepts: DMDA^using distributed arrays;
   Concepts: multicomponent
   Processors: n
T*/
/*F-----------------------------------------------------------------------

    We thank David E. Keyes for contributing the driven cavity discretization within this example code.

    This problem is modeled by the partial differential equation system

\begin{eqnarray}
        - \triangle U - \nabla_y \Omega & = & 0  \\
        - \triangle V + \nabla_x\Omega & = & 0  \\
        - \triangle \Omega + \nabla \cdot ([U*\Omega,V*\Omega]) - GR* \nabla_x T & = & 0  \\
        - \triangle T + PR* \nabla \cdot ([U*T,V*T]) & = & 0
\end{eqnarray}

    in the unit square, which is uniformly discretized in each of x and y in this simple encoding.

    No-slip, rigid-wall Dirichlet conditions are used for $ [U,V]$.
    Dirichlet conditions are used for Omega, based on the definition of
    vorticity: $ \Omega = - \nabla_y U + \nabla_x V$, where along each
    constant coordinate boundary, the tangential derivative is zero.
    Dirichlet conditions are used for T on the left and right walls,
    and insulation homogeneous Neumann conditions are used for T on
    the top and bottom walls.

    A finite difference approximation with the usual 5-point stencil
    is used to discretize the boundary value problem to obtain a
    nonlinear system of equations.  Upwinding is used for the divergence
    (convective) terms and central for the gradient (source) terms.

    The Jacobian can be either
      * formed via finite differencing using coloring (the default), or
      * applied matrix-free via the option -snes_mf
        (for larger grid problems this variant may not converge
        without a preconditioner due to ill-conditioning).

  ------------------------------------------------------------------------F*/

/*
   Include "petscdmda.h" so that we can use distributed arrays (DMDAs).
   Include "petscsnes.h" so that we can use SNES solvers.  Note that this
   file automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
     petscksp.h   - linear solvers
*/
#if defined(PETSC_APPLE_FRAMEWORK)
#import <PETSc/petscsnes.h>
#import <PETSc/petscdmda.h>
#else
#include <petscsnes.h>
#include <petscdm.h>
#include <petscdmda.h>
#endif
#include <iostream>
#include <utopia.hpp>      


/*
   User-defined routines and data structures
*/
typedef struct {
  PetscScalar u,v,omega,temp;
} Field;

PetscErrorCode FormFunctionLocal(DMDALocalInfo*,Field**,Field**,void*);

typedef struct {
  PetscReal   lidvelocity,prandtl,grashof;  /* physical parameters */
  PetscBool   draw_contours;                /* flag - 1 indicates drawing contours */
} AppCtx;

extern PetscErrorCode FormInitialGuess(AppCtx*,DM,Vec);
// extern PetscErrorCode NonlinearGS(SNES,Vec,Vec,void*);

extern PetscErrorCode ComputeMarker(DM da, Vec & Xguess);
extern PetscErrorCode Save_VTK_XML(DM da,Vec X,const char filename[]); 

int main(int argc,char **argv)
{
  AppCtx         user;                /* user-defined work context */
  PetscInt       mx,my,its;
  PetscErrorCode ierr;
  MPI_Comm       comm;
  SNES           snes;
  DM             da;
  Vec            x;

  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return(1);

  PetscFunctionBeginUser;
  comm = PETSC_COMM_WORLD;
  ierr = SNESCreate(comm,&snes);CHKERRQ(ierr);


  std::cout<<"inisde ---- \n"; 


  /*
      Create distributed array object to manage parallel grid and vectors
      for principal unknowns (x) and governing residuals (f)
  */
  ierr = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,4,4,PETSC_DECIDE,PETSC_DECIDE,4,1,0,0,&da);CHKERRQ(ierr);
  ierr = DMSetFromOptions(da);CHKERRQ(ierr);
  ierr = DMSetUp(da);CHKERRQ(ierr);
  ierr = SNESSetDM(snes,(DM)da);CHKERRQ(ierr);
  // ierr = SNESSetNGS(snes, NonlinearGS, (void*)&user);CHKERRQ(ierr);

  ierr = DMDAGetInfo(da,0,&mx,&my,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  /*
     Problem parameters (velocity of lid, prandtl, and grashof numbers)
  */
  user.lidvelocity = 1.0/(mx*my);
  user.prandtl     = 1.0;
  user.grashof     = 1.0;

  ierr = PetscOptionsGetReal(NULL,NULL,"-lidvelocity",&user.lidvelocity,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-prandtl",&user.prandtl,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-grashof",&user.grashof,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-contours",&user.draw_contours);CHKERRQ(ierr);

  ierr = DMDASetFieldName(da,0,"x_velocity");CHKERRQ(ierr);
  ierr = DMDASetFieldName(da,1,"y_velocity");CHKERRQ(ierr);
  ierr = DMDASetFieldName(da,2,"Omega");CHKERRQ(ierr);
  ierr = DMDASetFieldName(da,3,"temperature");CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create user context, set problem data, create vector data structures.
     Also, compute the initial guess.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create nonlinear solver context

     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMSetApplicationContext(da,&user);CHKERRQ(ierr);
  ierr = DMDASNESSetFunctionLocal(da,INSERT_VALUES,(PetscErrorCode (*)(DMDALocalInfo*,void*,void*,void*))FormFunctionLocal,&user);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  ierr = PetscPrintf(comm,"lid velocity = %g, prandtl # = %g, grashof # = %g\n",(double)user.lidvelocity,(double)user.prandtl,(double)user.grashof);CHKERRQ(ierr);


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve the nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
  ierr = FormInitialGuess(&user,da,x);CHKERRQ(ierr);

  // ierr = SNESSolve(snes,NULL,x);CHKERRQ(ierr);
  // ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);

  // SNESConvergedReason reason; 
  // ierr = SNESGetConvergedReason(snes, &reason); 
  // ierr = PetscPrintf(comm,"Number of SNES iterations = %D\n", its);CHKERRQ(ierr);
  // std::cout<<"reason: "<< reason << "  \n"; 


  Mat jac; 
  PetscInt dim; 
  VecGetSize(x, &dim); 

  ierr = MatCreate(PETSC_COMM_WORLD,&jac); CHKERRQ(ierr);
  ierr = MatSetSizes(jac, PETSC_DECIDE, PETSC_DECIDE, dim, dim); CHKERRQ(ierr);
  ierr = MatSetFromOptions(jac); CHKERRQ(ierr);
  ierr = MatSetUp(jac); CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,jac,jac,SNESComputeJacobianDefault,NULL); CHKERRQ(ierr);
  MatDestroy(&jac);

  Vec Marker; 
  std::cout<<"--- here.... \n"; 

  VecDuplicate(x, &Marker); 
  VecCopy(x, Marker); 

  ComputeMarker(da, Marker); 
  std::cout<<"--- out .... \n"; 

  
  {
    using namespace utopia; 
    DVectord x_u, Marker_u; 
    convert(x, x_u); 
    convert(Marker, Marker_u); 
    
    PETSCUtopiaNonlinearFunction<DSMatrixd, DVectord> fun(snes);

    auto linear_solver = std::make_shared<LUDecomposition<DSMatrixd, DVectord>>();
    // Newton<DSMatrixd, DVectord> newton(linear_solver); 
    // newton.verbose(true);  
    // newton.solve(fun, x_u); 

    // utopia::AffineSimilarity<utopia::DSMatrixd, utopia::DVectord> solver(linear_solver); 
    // AffineSimilarityAW<utopia::DSMatrixd, utopia::DVectord> solver(linear_solver); 

    utopia::PseudoContinuation<utopia::DSMatrixd, utopia::DVectord> solver(linear_solver); 


    //solver.set_scaling_matrix(utopia::local_identity(local_size(Mass_utopia).get(0), local_size(Mass_utopia).get(1))); 
    

      DSMatrixd Mass_utopia = identity(dim, dim); 
      std::cout<<"only in serial \n"; 
      std::cout<<"dim: "<< dim << "  \n"; 
      {
        Read<utopia::DVectord> read(Marker_u);
        Write<utopia::DSMatrixd> write(Mass_utopia);

        for(auto i=0; i < dim; i++)
        {
          if(Marker_u.get(i) == 1 ||Marker_u.get(i) == 2)
            Mass_utopia.set(i,i, 0);
        }
      }

    solver.set_mass_matrix(Mass_utopia); 
    solver.verbose(true);

    

    // solver.tau_option(3); 
    // solver.use_m(false); 
    // solver.set_m(-1); 
    solver.atol(1e-7); 
    solver.max_it(1000); 
    // solver.verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE); 
    solver.solve(fun, x_u); 



    convert(x_u, x); 
  }


  Save_VTK_XML(da, x, "output.vtk"); 

  /*
     Visualize solution
  */
  if (user.draw_contours) {
    ierr = VecView(x,PETSC_VIEWER_DRAW_WORLD);CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

/* ------------------------------------------------------------------- */

/*
   FormInitialGuess - Forms initial approximation.

   Input Parameters:
   user - user-defined application context
   X - vector

   Output Parameter:
   X - vector
*/
PetscErrorCode FormInitialGuess(AppCtx *user,DM da,Vec X)
{
  PetscInt       i,j,mx,xs,ys,xm,ym;
  PetscErrorCode ierr;
  PetscReal      grashof,dx;
  Field          **x;

  PetscFunctionBeginUser;
  grashof = user->grashof;

  ierr = DMDAGetInfo(da,0,&mx,0,0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  dx   = 1.0/(mx-1);

  /*
     Get local grid boundaries (for 2-dimensional DMDA):
       xs, ys   - starting grid indices (no ghost points)
       xm, ym   - widths of local grid (no ghost points)
  */
  ierr = DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);CHKERRQ(ierr);

  /*
     Get a pointer to vector data.
       - For default PETSc vectors, VecGetArray() returns a pointer to
         the data array.  Otherwise, the routine is implementation dependent.
       - You MUST call VecRestoreArray() when you no longer need access to
         the array.
  */
  ierr = DMDAVecGetArray(da,X,&x);CHKERRQ(ierr);

  /*
     Compute initial guess over the locally owned part of the grid
     Initial condition is motionless fluid and equilibrium temperature
  */
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      x[j][i].u     = 0.0;
      x[j][i].v     = 0.0;
      x[j][i].omega = 0.0;
      x[j][i].temp  = (grashof>0)*i*dx;
    }
  }

  /*
     Restore vector
  */
  ierr = DMDAVecRestoreArray(da,X,&x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FormFunctionLocal(DMDALocalInfo *info,Field **x,Field **f,void *ptr)
{
  AppCtx         *user = (AppCtx*)ptr;
  PetscErrorCode ierr;
  PetscInt       xints,xinte,yints,yinte,i,j;
  PetscReal      hx,hy,dhx,dhy,hxdhy,hydhx;
  PetscReal      grashof,prandtl,lid;
  PetscScalar    u,uxx,uyy,vx,vy,avx,avy,vxp,vxm,vyp,vym;

  PetscFunctionBeginUser;
  grashof = user->grashof;
  prandtl = user->prandtl;
  lid     = user->lidvelocity;

  /*
     Define mesh intervals ratios for uniform grid.

     Note: FD formulae below are normalized by multiplying through by
     local volume element (i.e. hx*hy) to obtain coefficients O(1) in two dimensions.


  */
  dhx   = (PetscReal)(info->mx-1);  dhy = (PetscReal)(info->my-1);
  hx    = 1.0/dhx;                   hy = 1.0/dhy;
  hxdhy = hx*dhy;                 hydhx = hy*dhx;

  xints = info->xs; xinte = info->xs+info->xm; yints = info->ys; yinte = info->ys+info->ym;

  /* Test whether we are on the bottom edge of the global array */
  if (yints == 0) {
    j     = 0;
    yints = yints + 1;
    /* bottom edge */
    for (i=info->xs; i<info->xs+info->xm; i++) {
      f[j][i].u     = x[j][i].u;
      f[j][i].v     = x[j][i].v;
      f[j][i].omega = x[j][i].omega + (x[j+1][i].u - x[j][i].u)*dhy;
      f[j][i].temp  = x[j][i].temp-x[j+1][i].temp;
    }
  }

  /* Test whether we are on the top edge of the global array */
  if (yinte == info->my) {
    j     = info->my - 1;
    yinte = yinte - 1;
    /* top edge */
    for (i=info->xs; i<info->xs+info->xm; i++) {
      f[j][i].u     = x[j][i].u - lid;
      f[j][i].v     = x[j][i].v;
      f[j][i].omega = x[j][i].omega + (x[j][i].u - x[j-1][i].u)*dhy;
      f[j][i].temp  = x[j][i].temp-x[j-1][i].temp;
    }
  }

  /* Test whether we are on the left edge of the global array */
  if (xints == 0) {
    i     = 0;
    xints = xints + 1;
    /* left edge */
    for (j=info->ys; j<info->ys+info->ym; j++) {
      f[j][i].u     = x[j][i].u;
      f[j][i].v     = x[j][i].v;
      f[j][i].omega = x[j][i].omega - (x[j][i+1].v - x[j][i].v)*dhx;
      f[j][i].temp  = x[j][i].temp;
    }
  }

  /* Test whether we are on the right edge of the global array */
  if (xinte == info->mx) {
    i     = info->mx - 1;
    xinte = xinte - 1;
    /* right edge */
    for (j=info->ys; j<info->ys+info->ym; j++) {
      f[j][i].u     = x[j][i].u;
      f[j][i].v     = x[j][i].v;
      f[j][i].omega = x[j][i].omega - (x[j][i].v - x[j][i-1].v)*dhx;
      f[j][i].temp  = x[j][i].temp - (PetscReal)(grashof>0);
    }
  }

  /* Compute over the interior points */
  for (j=yints; j<yinte; j++) {
    for (i=xints; i<xinte; i++) {

      /*
       convective coefficients for upwinding
      */
      vx  = x[j][i].u; avx = PetscAbsScalar(vx);
      vxp = .5*(vx+avx); vxm = .5*(vx-avx);
      vy  = x[j][i].v; avy = PetscAbsScalar(vy);
      vyp = .5*(vy+avy); vym = .5*(vy-avy);

      /* U velocity */
      u         = x[j][i].u;
      uxx       = (2.0*u - x[j][i-1].u - x[j][i+1].u)*hydhx;
      uyy       = (2.0*u - x[j-1][i].u - x[j+1][i].u)*hxdhy;
      f[j][i].u = uxx + uyy - .5*(x[j+1][i].omega-x[j-1][i].omega)*hx;

      /* V velocity */
      u         = x[j][i].v;
      uxx       = (2.0*u - x[j][i-1].v - x[j][i+1].v)*hydhx;
      uyy       = (2.0*u - x[j-1][i].v - x[j+1][i].v)*hxdhy;
      f[j][i].v = uxx + uyy + .5*(x[j][i+1].omega-x[j][i-1].omega)*hy;

      /* Omega */
      u             = x[j][i].omega;
      uxx           = (2.0*u - x[j][i-1].omega - x[j][i+1].omega)*hydhx;
      uyy           = (2.0*u - x[j-1][i].omega - x[j+1][i].omega)*hxdhy;
      f[j][i].omega = uxx + uyy + (vxp*(u - x[j][i-1].omega) + vxm*(x[j][i+1].omega - u))*hy +
                      (vyp*(u - x[j-1][i].omega) + vym*(x[j+1][i].omega - u))*hx -
                      .5*grashof*(x[j][i+1].temp - x[j][i-1].temp)*hy;

      /* Temperature */
      u            = x[j][i].temp;
      uxx          = (2.0*u - x[j][i-1].temp - x[j][i+1].temp)*hydhx;
      uyy          = (2.0*u - x[j-1][i].temp - x[j+1][i].temp)*hxdhy;
      f[j][i].temp =  uxx + uyy  + prandtl*((vxp*(u - x[j][i-1].temp) + vxm*(x[j][i+1].temp - u))*hy +
                                            (vyp*(u - x[j-1][i].temp) + vym*(x[j+1][i].temp - u))*hx);
    }
  }

  /*
     Flop count (multiply-adds are counted as 2 operations)
  */
  ierr = PetscLogFlops(84.0*info->ym*info->xm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}



PetscErrorCode ComputeMarker(DM da, Vec & Xguess)
/* ------------------------------------------------------------------- */
{
  AppCtx         *user;
  PetscInt       i,j,is,js,im,jm;
  PetscErrorCode ierr;
  Field          **x_new;

  /* Get the fine grid */
  ierr   = DMGetApplicationContext(da,&user);CHKERRQ(ierr);
  ierr   = DMDAGetCorners(da,&is,&js,NULL,&im,&jm,NULL);CHKERRQ(ierr);
  ierr   = DMDAVecGetArray(da,Xguess,(void**)&x_new);CHKERRQ(ierr);

  /* Compute initial guess */
  for (j=js; j<js+jm; j++) {
    for (i=is; i<is+im; i++) {
      x_new[j][i].u = 1.0; 
      x_new[j][i].v = 2.0;
      x_new[j][i].omega = 3.0; 
      x_new[j][i].temp = 4.0; 
    }
  }

  /* Restore x to Xguess */
  ierr = DMDAVecRestoreArray(da,Xguess,(void**)&x_new);CHKERRQ(ierr);

  return 0;
}



PetscErrorCode Save_VTK_XML(DM da,Vec X,const char filename[])
{
  PetscErrorCode ierr;
  PetscViewer       viewer;
  // ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&view);CHKERRQ(ierr);
  // ierr = VecView(X,view);CHKERRQ(ierr);
  // ierr = PetscViewerDestroy(&view);CHKERRQ(ierr);


  PetscViewerASCIIOpen(PETSC_COMM_WORLD, "testing.vtk", &viewer);
  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK);
  DMDASetUniformCoordinates(da, 0.0, 1, 0.0, 1, 0.0, 0.0);
  DMView(da, viewer);
  VecView(X, viewer);
  PetscViewerDestroy(&viewer);



  PetscFunctionReturn(0);
}
