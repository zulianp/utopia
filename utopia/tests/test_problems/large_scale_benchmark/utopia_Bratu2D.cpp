#include "utopia_Bratu2D.hpp"
#include "utopia_TestFunctions.hpp"

#ifdef WITH_PETSC

namespace utopia {

    PetscErrorCode Bratu2DFormFunctionLocal(DMDALocalInfo *info,
                                            PetscScalar **x,
                                            PetscScalar **f,
                                            ParamsBratu2D *user) {
        PetscErrorCode ierr;
        PetscInt i, j;
        PetscReal lambda, hx, hy, hxdhy, hydhx;
        PetscScalar u, ue, uw, un, us, uxx, uyy, mms_solution;
        DMDACoor2d c;

        PetscFunctionBeginUser;
        lambda = user->lambda;
        hx = 1.0 / static_cast<PetscReal>(info->mx - 1);
        hy = 1.0 / static_cast<PetscReal>(info->my - 1);
        hxdhy = hx / hy;
        hydhx = hy / hx;
        /*
           Compute function over the locally owned part of the grid
        */
        for (j = info->ys; j < info->ys + info->ym; j++) {
            for (i = info->xs; i < info->xs + info->xm; i++) {
                if (i == 0 || j == 0 || i == info->mx - 1 || j == info->my - 1) {
                    c.x = i * hx;
                    c.y = j * hy;
                    ierr = user->mms_solution(user, &c, &mms_solution);
                    CHKERRQ(ierr);
                    f[j][i] = 2.0 * (hydhx + hxdhy) * (x[j][i] - mms_solution);
                } else {
                    u = x[j][i];
                    uw = x[j][i - 1];
                    ue = x[j][i + 1];
                    un = x[j - 1][i];
                    us = x[j + 1][i];

                    /* Enforce boundary conditions at neighboring points -- setting these values causes the Jacobian to
                     * be symmetric. */
                    if (i - 1 == 0) {
                        c.x = (i - 1) * hx;
                        c.y = j * hy;
                        ierr = user->mms_solution(user, &c, &uw);
                        CHKERRQ(ierr);
                    }
                    if (i + 1 == info->mx - 1) {
                        c.x = (i + 1) * hx;
                        c.y = j * hy;
                        ierr = user->mms_solution(user, &c, &ue);
                        CHKERRQ(ierr);
                    }
                    if (j - 1 == 0) {
                        c.x = i * hx;
                        c.y = (j - 1) * hy;
                        ierr = user->mms_solution(user, &c, &un);
                        CHKERRQ(ierr);
                    }
                    if (j + 1 == info->my - 1) {
                        c.x = i * hx;
                        c.y = (j + 1) * hy;
                        ierr = user->mms_solution(user, &c, &us);
                        CHKERRQ(ierr);
                    }

                    uxx = (2.0 * u - uw - ue) * hydhx;
                    uyy = (2.0 * u - un - us) * hxdhy;
                    c.x = i * hx;
                    c.y = j * hy;
                    f[j][i] = uxx + uyy - hx * hy * (lambda * std::exp(u));
                }
            }
        }
        ierr = PetscLogFlops(11.0 * info->ym * info->xm);
        CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }

    PetscErrorCode Bratu2DFormJacobianLocal(DMDALocalInfo *info,
                                            PetscScalar **x,
                                            Mat jac,
                                            Mat jacpre,
                                            ParamsBratu2D *user) {
        PetscErrorCode ierr;
        PetscInt i, j, k;
        MatStencil col[5], row;
        PetscScalar lambda, v[5], hx, hy, hxdhy, hydhx, sc;
        DM coordDA;
        Vec coordinates;
        DMDACoor2d **coords;

        PetscFunctionBeginUser;
        lambda = user->lambda;
        /* Extract coordinates */
        ierr = DMGetCoordinateDM(info->da, &coordDA);
        CHKERRQ(ierr);
        ierr = DMGetCoordinates(info->da, &coordinates);
        CHKERRQ(ierr);
        ierr = DMDAVecGetArray(coordDA, coordinates, &coords);
        CHKERRQ(ierr);
        hx = info->xm > 1
                 ? PetscRealPart(coords[info->ys][info->xs + 1].x) - PetscRealPart(coords[info->ys][info->xs].x)
                 : 1.0;
        hy = info->ym > 1
                 ? PetscRealPart(coords[info->ys + 1][info->xs].y) - PetscRealPart(coords[info->ys][info->xs].y)
                 : 1.0;
        ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);
        CHKERRQ(ierr);
        hxdhy = hx / hy;
        hydhx = hy / hx;
        sc = hx * hy * lambda;

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
        for (j = info->ys; j < info->ys + info->ym; j++) {
            for (i = info->xs; i < info->xs + info->xm; i++) {
                row.j = j;
                row.i = i;
                /* boundary points */
                if (i == 0 || j == 0 || i == info->mx - 1 || j == info->my - 1) {
                    v[0] = 2.0 * (hydhx + hxdhy);
                    ierr = MatSetValuesStencil(jacpre, 1, &row, 1, &row, v, INSERT_VALUES);
                    CHKERRQ(ierr);
                } else {
                    k = 0;
                    /* interior grid points */
                    if (j - 1 != 0) {
                        v[k] = -hxdhy;
                        col[k].j = j - 1;
                        col[k].i = i;
                        k++;
                    }
                    if (i - 1 != 0) {
                        v[k] = -hydhx;
                        col[k].j = j;
                        col[k].i = i - 1;
                        k++;
                    }

                    v[k] = 2.0 * (hydhx + hxdhy) - sc * std::exp(x[j][i]);
                    col[k].j = row.j;
                    col[k].i = row.i;
                    k++;

                    if (i + 1 != info->mx - 1) {
                        v[k] = -hydhx;
                        col[k].j = j;
                        col[k].i = i + 1;
                        k++;
                    }
                    if (j + 1 != info->mx - 1) {
                        v[k] = -hxdhy;
                        col[k].j = j + 1;
                        col[k].i = i;
                        k++;
                    }
                    ierr = MatSetValuesStencil(jacpre, 1, &row, k, col, v, INSERT_VALUES);
                    CHKERRQ(ierr);
                }
            }
        }

        /*
           Assemble matrix, using the 2-step process:
             MatAssemblyBegin(), MatAssemblyEnd().
        */
        ierr = MatAssemblyBegin(jacpre, MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
        ierr = MatAssemblyEnd(jacpre, MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
        /*
           Tell the matrix we will never add a new nonzero location to the
           matrix. If we do, it will generate an error.
        */
        ierr = MatSetOption(jac, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
        CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }

    PetscErrorCode Bratu2DFormObjectiveLocal(DMDALocalInfo *info,
                                             PetscScalar **x,
                                             PetscReal *obj,
                                             ParamsBratu2D *user) {
        PetscErrorCode ierr;
        PetscInt i, j;
        PetscReal lambda, hx, hy, hxdhy, hydhx, sc, lobj = 0;
        PetscScalar u, ue, uw, un, us, uxux, uyuy;
        MPI_Comm comm;

        PetscFunctionBeginUser;
        *obj = 0;
        ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(info->da), &comm);
        CHKERRQ(ierr);
        lambda = user->lambda;
        hx = 1.0 / static_cast<PetscReal>(info->mx - 1);
        hy = 1.0 / static_cast<PetscReal>(info->my - 1);
        sc = hx * hy * lambda;
        hxdhy = hx / hy;
        hydhx = hy / hx;
        /*
           Compute function over the locally owned part of the grid
        */
        for (j = info->ys; j < info->ys + info->ym; j++) {
            for (i = info->xs; i < info->xs + info->xm; i++) {
                if (i == 0 || j == 0 || i == info->mx - 1 || j == info->my - 1) {
                    lobj += PetscRealPart((hydhx + hxdhy) * x[j][i] * x[j][i]);
                    // lobj += PetscRealPart((hydhx + hxdhy)*x[j][i]*x[j][i]) - sc*std::exp(u);
                } else {
                    u = x[j][i];
                    uw = x[j][i - 1];
                    ue = x[j][i + 1];
                    un = x[j - 1][i];
                    us = x[j + 1][i];

                    if (i - 1 == 0) {
                        uw = 0.;
                    }
                    if (i + 1 == info->mx - 1) {
                        ue = 0.;
                    }
                    if (j - 1 == 0) {
                        un = 0.;
                    }
                    if (j + 1 == info->my - 1) {
                        us = 0.;
                    }

                    /* F[u] = 1/2\int_{\omega}\nabla^2u(x)*u(x)*dx */
                    uxux = u * (2. * u - ue - uw) * hydhx;
                    uyuy = u * (2. * u - un - us) * hxdhy;
                    lobj += PetscRealPart(0.5 * (uxux + uyuy) - sc * std::exp(u));
                }
            }
        }
        ierr = PetscLogFlops(12.0 * info->ym * info->xm);
        CHKERRQ(ierr);
        ierr = MPI_Allreduce(&lobj, obj, 1, MPIU_REAL, MPIU_SUM, comm);
        CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }

    PetscErrorCode Bratu2DMMSSolution(ParamsBratu2D * /*user*/, const DMDACoor2d *c, PetscScalar *u) {
        // we need to return 0 for bc
        PetscReal x = PetscRealPart(c->x), y = PetscRealPart(c->y);
        u[0] = 0.0;
        return 0;
    }

    PetscErrorCode Bratu2DFormBCData(DM da, ParamsBratu2D * /*user*/, Vec BC_flag, Vec BC_value) {
        PetscInt i, j, Mx, My, xs, ys, xm, ym;
        PetscErrorCode ierr;
        PetscScalar **BC_value_w;
        PetscScalar **BC_flag_w;
        DM coordDA;
        Vec coordinates;
        DMDACoor2d **coords;

        PetscFunctionBeginUser;
        ierr = DMDAGetInfo(da,
                           PETSC_IGNORE,
                           &Mx,
                           &My,
                           PETSC_IGNORE,
                           PETSC_IGNORE,
                           PETSC_IGNORE,
                           PETSC_IGNORE,
                           PETSC_IGNORE,
                           PETSC_IGNORE,
                           PETSC_IGNORE,
                           PETSC_IGNORE,
                           PETSC_IGNORE,
                           PETSC_IGNORE);
        CHKERRQ(ierr);

        ierr = DMDAVecGetArray(da, BC_value, &BC_value_w);
        CHKERRQ(ierr);
        ierr = DMDAVecGetArray(da, BC_flag, &BC_flag_w);
        CHKERRQ(ierr);

        ierr = DMDAGetCorners(da, &xs, &ys, nullptr, &xm, &ym, nullptr);
        CHKERRQ(ierr);

        ierr = DMGetCoordinateDM(da, &coordDA);
        CHKERRQ(ierr);
        ierr = DMGetCoordinates(da, &coordinates);
        CHKERRQ(ierr);
        ierr = DMDAVecGetArray(coordDA, coordinates, &coords);
        CHKERRQ(ierr);

        /*
           Compute initial guess over the locally owned part of the grid
        */
        for (j = ys; j < ys + ym; j++) {
            for (i = xs; i < xs + xm; i++) {
                if (i == 0 || j == 0 || i == Mx - 1 || j == My - 1) {
                    /* boundary conditions are all zero Dirichlet */
                    BC_flag_w[j][i] = 1.0;
                    BC_value_w[j][i] = 0.0;

                } else {
                    BC_flag_w[j][i] = 0.0;
                    BC_value_w[j][i] = 0.0;
                }
            }
        }

        /*
           Restore vector
        */
        ierr = DMDAVecRestoreArray(da, BC_value, &BC_value_w);
        CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(da, BC_flag, &BC_flag_w);
        CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);
        CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }

    PetscErrorCode Bratu2DFormInitialGuess(DM da, ParamsBratu2D *user, Vec X) {
        PetscInt i, j, Mx, My, xs, ys, xm, ym;
        PetscErrorCode ierr;
        PetscReal lambda, temp1, temp, hx, hy;
        PetscScalar **x;

        PetscFunctionBeginUser;
        ierr = DMDAGetInfo(da,
                           PETSC_IGNORE,
                           &Mx,
                           &My,
                           PETSC_IGNORE,
                           PETSC_IGNORE,
                           PETSC_IGNORE,
                           PETSC_IGNORE,
                           PETSC_IGNORE,
                           PETSC_IGNORE,
                           PETSC_IGNORE,
                           PETSC_IGNORE,
                           PETSC_IGNORE,
                           PETSC_IGNORE);
        CHKERRQ(ierr);

        lambda = user->lambda;
        hx = 1.0 / static_cast<PetscReal>(Mx - 1);
        hy = 1.0 / static_cast<PetscReal>(My - 1);
        temp1 = lambda / (lambda + 1.0);

        /*
           Get a pointer to vector data.
             - For default PETSc vectors, VecGetArray() returns a pointer to
               the data array.  Otherwise, the routine is implementation dependent.
             - You MUST call VecRestoreArray() when you no longer need access to
               the array.
        */
        ierr = DMDAVecGetArray(da, X, &x);
        CHKERRQ(ierr);

        /*
           Get local grid boundaries (for 2-dimensional DMDA):
             xs, ys   - starting grid indices (no ghost points)
             xm, ym   - widths of local grid (no ghost points)

        */
        ierr = DMDAGetCorners(da, &xs, &ys, nullptr, &xm, &ym, nullptr);
        CHKERRQ(ierr);

        /*
           Compute initial guess over the locally owned part of the grid
        */
        for (j = ys; j < ys + ym; j++) {
            temp = static_cast<PetscReal>(PetscMin(j, My - j - 1)) * hy;
            for (i = xs; i < xs + xm; i++) {
                if (i == 0 || j == 0 || i == Mx - 1 || j == My - 1) {
                    /* boundary conditions are all zero Dirichlet */
                    x[j][i] = 0.0;
                } else {
                    x[j][i] = temp1 * std::sqrt(PetscMin((PetscReal)(PetscMin(i, Mx - i - 1)) * hx, temp));
                }
            }
        }

        /*
           Restore vector
        */
        ierr = DMDAVecRestoreArray(da, X, &x);
        CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }

    // explicit
    template class Bratu2D<PetscMatrix, PetscVector, PETSC>;
}  // namespace utopia

#endif  // WITH_PETSC
