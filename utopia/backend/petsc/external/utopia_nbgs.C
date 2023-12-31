#include <math.h>
#include <utopia_nbgs.h>
#include <petsc/private/kernels/blockinvert.h>
#include <stdio.h>
#include <petscmat.h>
//#include <aij.h>

#include "aij.h"
#include <petscsys.h>

/*

 Placeholder for box constraints

 */
#undef __FUNCT__
#define __FUNCT__ "constrain_box"
PetscErrorCode  constrain_box(BLOCKSOLVER * blocksolver, MatScalar * diag, PetscScalar * t, PetscScalar * x,
                              PetscScalar * lb, PetscScalar * ub, int row,
                              int bs, PetscBool * constrained){
    PetscFunctionBegin;

    for(int i = 0; i<bs;i++){
        constrained[i] = PETSC_FALSE;
        if (lb[row + i]> ub[row + i]){
            PetscPrintf(PETSC_COMM_WORLD, "Constraints not well-defined: Lower boundary bigger than upper boundary \n");
            PetscFunctionReturn(0);
        }
        if (x[row + i] < lb[row + i] + NBGS_TOL){
            x[row + i] = lb[row + i];
            constrained[i] = PETSC_TRUE;
        }
        if (x[row + i] > ub[row + i] + NBGS_TOL){
            x[row + i] = ub[row + i];
            constrained[i] = PETSC_TRUE;
        }
    };
    PetscFunctionReturn(0);
};

#undef __FUNCT__
#define __FUNCT__ "constrain_blocknd_ub"
PetscErrorCode  constrain_blocknd_ub(BLOCKSOLVER * blocksolver, MatScalar * diag, PetscScalar * t, PetscScalar * x,
                                  PetscScalar * lb, PetscScalar * ub, int row,
                                  int bs, PetscBool * constrained){
    PetscFunctionBegin;
    PetscErrorCode ierr;

//    MatScalar * diag_constrained;
    PetscBool allowzeropivot=PETSC_TRUE,zeropivotdetected=PETSC_FALSE;

//    PetscMalloc1(bs*bs, &blocksolver->diag_constrained);

    for(int i = 0; i<bs;i++){
        constrained[i] = PETSC_FALSE;
        for(int j = 0; j<bs;j++){
            blocksolver->diag_constrained[i*bs+j] = diag[i*bs+j];
        }
    }

    for(int i = 0; i<bs;i++){
        if (x[row + i] > ub[row + i]+ NBGS_TOL){
            for (int j=0;j<bs;j++){
                blocksolver->diag_constrained[j*bs+i] = 0;
            }
            constrained[i] = PETSC_TRUE;
            blocksolver->diag_constrained[i*bs+i] = 1;
            t[row+i]=ub[row+i];

            if (bs == 2){
                ierr = PetscKernel_A_gets_inverse_A_2(blocksolver->diag_constrained,0,allowzeropivot,&zeropivotdetected);

                x[row] = t[row]*blocksolver->diag_constrained[0] + t[row+1]*blocksolver->diag_constrained[2];
                x[row+1] = t[row]*blocksolver->diag_constrained[1] + t[row+1]*blocksolver->diag_constrained[3] ;
            }
            else if (bs == 3){
                ierr = PetscKernel_A_gets_inverse_A_3(blocksolver->diag_constrained,0,allowzeropivot,&zeropivotdetected);

                x[row] = t[row]*blocksolver->diag_constrained[0] + t[row+1]*blocksolver->diag_constrained[3] + t[row+2]*blocksolver->diag_constrained[6];
                x[row+1] = t[row]*blocksolver->diag_constrained[1] + t[row+1]*blocksolver->diag_constrained[4] + t[row+2]*blocksolver->diag_constrained[7];
                x[row+2] = t[row]*blocksolver->diag_constrained[2] + t[row+1]*blocksolver->diag_constrained[5] + t[row+2]*blocksolver->diag_constrained[8];

            }else if (bs == 1){
                x[row] = t[row]/blocksolver->diag_constrained[0];
            }else{
                PetscPrintf(PETSC_COMM_WORLD, "NBGS: Blocksize %i is not supported. Exiting. \n", bs);
                return(0);
            }
        };
    };


    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "constrain_blocknd_lb"
PetscErrorCode  constrain_blocknd_lb(BLOCKSOLVER * blocksolver, MatScalar * diag, PetscScalar * t, PetscScalar * x,
                                  PetscScalar * lb, PetscScalar * ub, int row,
                                  int bs, PetscBool * constrained){
    PetscFunctionBegin;

    PetscErrorCode ierr;

//    MatScalar * diag_constrained;
    PetscBool allowzeropivot=PETSC_TRUE,zeropivotdetected=PETSC_FALSE;

//    PetscMalloc1(bs*bs, &blocksolver->diag_constrained);

    for(int i = 0; i<bs;i++){
        constrained[i] = PETSC_FALSE;
        for(int j = 0; j<bs;j++){
            blocksolver->diag_constrained[i*bs+j] = diag[i*bs+j];
        }
    }

    for(int i = 0; i<bs;i++){
        if (x[row + i] < lb[row + i]+NBGS_TOL){
            for (int j=0;j<bs;j++){
                blocksolver->diag_constrained[j*bs+i] = 0;
            }
            constrained[i] = PETSC_TRUE;
            blocksolver->diag_constrained[i*bs+i] = 1;
            t[row+i]=lb[row+i];

            if (bs == 2){
                ierr = PetscKernel_A_gets_inverse_A_2(blocksolver->diag_constrained,0,allowzeropivot,&zeropivotdetected);

                x[row] = t[row]*blocksolver->diag_constrained[0] + t[row+1]*blocksolver->diag_constrained[2];
                x[row+1] = t[row]*blocksolver->diag_constrained[1] + t[row+1]*blocksolver->diag_constrained[3] ;
            }
            else if (bs == 3){
                ierr = PetscKernel_A_gets_inverse_A_3(blocksolver->diag_constrained,0,allowzeropivot,&zeropivotdetected);

                x[row] = t[row]*blocksolver->diag_constrained[0] + t[row+1]*blocksolver->diag_constrained[3] + t[row+2]*blocksolver->diag_constrained[6];
                x[row+1] = t[row]*blocksolver->diag_constrained[1] + t[row+1]*blocksolver->diag_constrained[4] + t[row+2]*blocksolver->diag_constrained[7];
                x[row+2] = t[row]*blocksolver->diag_constrained[2] + t[row+1]*blocksolver->diag_constrained[5] + t[row+2]*blocksolver->diag_constrained[8];

            }else if (bs == 1){
                x[row] = t[row]/blocksolver->diag_constrained[0];
            } else {
                x[row] = t[row]/diag[0];
                PetscPrintf(PETSC_COMM_WORLD, "NBGS: Blocksize %i is not supported. Exiting. \n", bs);
                return(0);
            }
        };
    };
//    PetscFree(diag_constrained);
    PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "constrain_blocknd"
PetscErrorCode  constrain_blocknd(BLOCKSOLVER * blocksolver, MatScalar * diag, PetscScalar * t, PetscScalar * x,
                                  PetscScalar * lb, PetscScalar * ub, int row,
                                  int bs, PetscBool * constrained){
    PetscFunctionBegin;

    PetscErrorCode ierr;

//    MatScalar * diag_constrained;
    PetscBool allowzeropivot=PETSC_TRUE,zeropivotdetected=PETSC_FALSE;

//    PetscMalloc1(bs*bs, &blocksolver->diag_constrained);

    for(int i = 0; i<bs;i++){
        constrained[i] = PETSC_FALSE;
        for(int j = 0; j<bs;j++){
            blocksolver->diag_constrained[i*bs+j] = diag[i*bs+j];
        }
        if (lb[row + i]> ub[row + i]){
            PetscPrintf(PETSC_COMM_WORLD, "Constraints not well-defined: Lower boundary bigger than upper boundary \n");
            PetscFunctionReturn(0);
        }
    }

    for(int i = 0; i<bs;i++){
        if (x[row + i] < lb[row + i]+NBGS_TOL){
            for (int j=0;j<bs;j++){
                blocksolver->diag_constrained[j*bs+i] = 0;
            }
            constrained[i] = PETSC_TRUE;
            blocksolver->diag_constrained[i*bs+i] = 1;
            t[row+i]=lb[row+i];

            if (bs == 2){
                ierr = PetscKernel_A_gets_inverse_A_2(blocksolver->diag_constrained,0,allowzeropivot,&zeropivotdetected);

                x[row] = t[row]*blocksolver->diag_constrained[0] + t[row+1]*blocksolver->diag_constrained[2];
                x[row+1] = t[row]*blocksolver->diag_constrained[1] + t[row+1]*blocksolver->diag_constrained[3] ;
            }
            else if (bs == 3){
                ierr = PetscKernel_A_gets_inverse_A_3(blocksolver->diag_constrained,0,allowzeropivot,&zeropivotdetected);

                x[row] = t[row]*blocksolver->diag_constrained[0] + t[row+1]*blocksolver->diag_constrained[3] + t[row+2]*blocksolver->diag_constrained[6];
                x[row+1] = t[row]*blocksolver->diag_constrained[1] + t[row+1]*blocksolver->diag_constrained[4] + t[row+2]*blocksolver->diag_constrained[7];
                x[row+2] = t[row]*blocksolver->diag_constrained[2] + t[row+1]*blocksolver->diag_constrained[5] + t[row+2]*blocksolver->diag_constrained[8];

            }else if (bs == 1){
                x[row] = t[row]/blocksolver->diag_constrained[0];
            } else {
                x[row] = t[row]/diag[0];
                PetscPrintf(PETSC_COMM_WORLD, "NBGS: Blocksize %i is not supported. Exiting. \n", bs);
                return(0);
            }
        };
    };

    for(int i = 0; i<bs;i++){
        if (x[row + i] > ub[row + i]+ NBGS_TOL){
            for (int j=0;j<bs;j++){
                blocksolver->diag_constrained[j*bs+i] = 0;
            }
            constrained[i] = PETSC_TRUE;
            blocksolver->diag_constrained[i*bs+i] = 1;
            t[row+i]=ub[row+i];

            if (bs == 2){
                ierr = PetscKernel_A_gets_inverse_A_2(blocksolver->diag_constrained,0,allowzeropivot,&zeropivotdetected);

                x[row] = t[row]*blocksolver->diag_constrained[0] + t[row+1]*blocksolver->diag_constrained[2];
                x[row+1] = t[row]*blocksolver->diag_constrained[1] + t[row+1]*blocksolver->diag_constrained[3] ;
            }
            else if (bs == 3){
                ierr = PetscKernel_A_gets_inverse_A_3(blocksolver->diag_constrained,0,allowzeropivot,&zeropivotdetected);

                x[row] = t[row]*blocksolver->diag_constrained[0] + t[row+1]*blocksolver->diag_constrained[3] + t[row+2]*blocksolver->diag_constrained[6];
                x[row+1] = t[row]*blocksolver->diag_constrained[1] + t[row+1]*blocksolver->diag_constrained[4] + t[row+2]*blocksolver->diag_constrained[7];
                x[row+2] = t[row]*blocksolver->diag_constrained[2] + t[row+1]*blocksolver->diag_constrained[5] + t[row+2]*blocksolver->diag_constrained[8];

            }else{
                PetscPrintf(PETSC_COMM_WORLD, "NBGS: Blocksize %i is not supported. Exiting. \n", bs);
                return(0);
            }
        };
    };

//    PetscFree(diag_constrained);
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "local_inode_nbgs"
// Most of this is stolen from the petsc-source. Only works for aij-matrices.
PetscErrorCode local_inode_nbgs(NBGS_CTX * nbgs, Mat A,Vec bb, Vec xx, Vec lblb, Vec ubub, int bs, PetscReal omega){

    Mat_SeqAIJ        *a = (Mat_SeqAIJ*)A->data;
    PetscScalar       sum1 = 0.0,sum2 = 0.0,sum3 = 0.0,sum4 = 0.0,sum5 = 0.0,tmp0,tmp1;
    MatScalar         *ibdiag,*bdiag,work[25],*t;
    PetscScalar       *x, *lb, *ub;
    const MatScalar   *v = a->a,*v1 = NULL,*v2 = NULL,*v3 = NULL,*v4 = NULL,*v5 = NULL;
    const PetscScalar *xb, *b;
    PetscReal         zeropivot = 1.0e-15, shift = 0.0;
    PetscErrorCode    ierr;
    PetscInt          n,m = a->inode.node_count,cnt = 0,i,j,row,i1,i2,nodes;
    PetscInt          sz,k,ipvt[5];
    const PetscInt    *sizes = a->inode.size,*idx,*diag = a->diag,*ii = a->i;
    PetscBool         allowzeropivot,zeropivotdetected=PETSC_FALSE;
    PetscBool         *constrained_nodal_coords;

    //    Vec loc_active_indices;

    PetscFunctionBegin;

    ierr = PetscMalloc(bs, &constrained_nodal_coords);CHKERRQ(ierr);

    ierr = MatGetLocalSize(A, &m, &n); CHKERRQ(ierr);
    nodes = m/bs;
    // we need some security checks here for non-fitting inode-setups ...
    if (sizes != 0 && bs >1){
        for (i=0; i<nodes;i++){
            if (sizes[i] != bs ){
                PetscPrintf(PETSC_COMM_WORLD,"Warning (NBGS): Inode sizes don't match block sizes (node=%i, size[i]=%i, bs=%i) \n", i, sizes[i], bs);
            }
        }
    }

    if (nbgs->constrained && nbgs->active_coordinates){
        //        ierr = VecDuplicate(xx, &loc_active_indices); CHKERRQ(ierr);
        ierr = VecGetLocalVector(nbgs->active_coordinates, nbgs->loc_active_indices); CHKERRQ(ierr);
    }

    if (!a->inode.ibdiagvalid) {
        if (!a->inode.ibdiag) {
            /* calculate space needed for diagonal blocks */
            for (i=0; i<nodes; i++) {
                cnt += bs*bs; // we essentially overwrite the original inode sizes with blocksize
            }
            a->inode.bdiagsize = cnt;
            //            ierr = PetscMalloc3(cnt,&a->inode.ibdiag,cnt,&a->inode.bdiag,A->rmap->n,&a->inode.ssor_work);CHKERRQ(ierr);
            ierr = PetscMalloc1(cnt,&a->inode.ibdiag); CHKERRQ(ierr);
            ierr = PetscMalloc1(cnt,&a->inode.bdiag); CHKERRQ(ierr);
            ierr = PetscMalloc1(A->rmap->n,&a->inode.ssor_work); CHKERRQ(ierr);
        }

        /* copy over the diagonal blocks and invert them */
        ibdiag = a->inode.ibdiag;
        bdiag  = a->inode.bdiag;
        cnt    = 0;
        for (i=0, row = 0; i<nodes; i++) {
            for (j=0; j<bs; j++) {
                for (k=0; k<bs; k++) {
                    bdiag[cnt+k*bs+j] = v[diag[row+j] - j + k];
                }
            }
            ierr = PetscMemcpy(ibdiag+cnt,bdiag+cnt,bs*bs*sizeof(MatScalar));CHKERRQ(ierr);
            allowzeropivot = PETSC_TRUE;
            switch (bs) {
                case 1:
                    /* Create matrix data structure */
                    if (PetscAbsScalar(ibdiag[cnt]) < zeropivot) PetscPrintf(PETSC_COMM_WORLD,"Fixme (local_inode_nbgs): Zero pivot on row %D",row);
                    ibdiag[cnt] = 1.0/ibdiag[cnt];
                    break;
                case 2:
                    ierr = PetscKernel_A_gets_inverse_A_2(ibdiag+cnt,shift,allowzeropivot,&zeropivotdetected);CHKERRQ(ierr);
                    break;
                case 3:
                    ierr = PetscKernel_A_gets_inverse_A_3(ibdiag+cnt,shift,allowzeropivot,&zeropivotdetected);CHKERRQ(ierr);
                    break;
                case 4:
                    ierr = PetscKernel_A_gets_inverse_A_4(ibdiag+cnt,shift,allowzeropivot,&zeropivotdetected);CHKERRQ(ierr);
                    break;
                case 5:
                    ierr = PetscKernel_A_gets_inverse_A_5(ibdiag+cnt,ipvt,work,shift,allowzeropivot,&zeropivotdetected);CHKERRQ(ierr);
                    break;
                default:
                    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"Inode size %D not supported",sizes[i]);
            }
            cnt += bs*bs;
            row += bs;
        }
        a->inode.ibdiagvalid = PETSC_TRUE;
    }
    ibdiag = a->inode.ibdiag;
    bdiag  = a->inode.bdiag;
    t      = a->inode.ssor_work;

    ierr = VecGetArray(xx,&x);CHKERRQ(ierr);
    ierr = VecGetArrayRead(bb,&b);CHKERRQ(ierr);
    if (nbgs->constrained) {
        ierr = VecGetArray(lblb,&lb);CHKERRQ(ierr);
        ierr = VecGetArray(ubub,&ub);CHKERRQ(ierr);
    }

    /* We count flops by assuming the upper triangular and lower triangular parts have the same number of nonzeros */
    // SOR_ZERO_INITIAL_GUESS)
    // SOR_FORWARD_SWEEP || flag & SOR_LOCAL_FORWARD_SWEEP
    for (i=0, row=0; i<nodes; i++) {
        sz  = diag[row] - ii[row];
        v1  = a->a + ii[row];
        idx = a->j + ii[row];

        /* see comments for MatMult_SeqAIJ_Inode() for how this is coded */
        switch (bs){
            case 1:

                sum1 = b[row];
                for (n = 0; n<sz-1; n+=2) {
                    i1    = idx[0];
                    i2    = idx[1];
                    idx  += 2;
                    tmp0  = x[i1];
                    tmp1  = x[i2];
                    sum1 -= v1[0] * tmp0 + v1[1] * tmp1; v1 += 2;
                }

                if (n == sz-1) {
                    tmp0  = x[*idx];
                    sum1 -= *v1 * tmp0;
                }
                t[row]   = sum1;
                x[row++] = sum1*(*ibdiag++);

                if (nbgs->constrained){
                    nbgs->constrain(&nbgs->blocksolver, bdiag, t, x, lb, ub, (i*bs), bs, constrained_nodal_coords);
                    for (int j=0;j<bs;j++){
                        if(constrained_nodal_coords[j]){
                            ierr = VecSetValue(nbgs->loc_active_indices,(i*bs), 1, INSERT_VALUES); CHKERRQ(ierr);

//                            ierr = VecAssemblyBegin(nbgs->loc_active_indices);CHKERRQ(ierr);
//                            ierr = VecAssemblyEnd(nbgs->loc_active_indices);CHKERRQ(ierr);
                        }
                    }
                }

                break;
            case 2:
                v2   = a->a + ii[row+1];
                sum1 = b[row];
                sum2 = b[row+1];
                for (n = 0; n<sz-1; n+=2) {
                    i1    = idx[0];
                    i2    = idx[1];
                    idx  += 2;
                    tmp0  = x[i1];
                    tmp1  = x[i2];
                    sum1 -= v1[0] * tmp0 + v1[1] * tmp1; v1 += 2;
                    sum2 -= v2[0] * tmp0 + v2[1] * tmp1; v2 += 2;
                }

                if (n == sz-1) {
                    tmp0  = x[*idx];
                    sum1 -= v1[0] * tmp0;
                    sum2 -= v2[0] * tmp0;
                }
                t[row]   = sum1;
                t[row+1] = sum2;
                x[row++] = sum1*ibdiag[0] + sum2*ibdiag[2];
                x[row++] = sum1*ibdiag[1] + sum2*ibdiag[3];


                if (nbgs->constrained){
                    nbgs->constrain(&nbgs->blocksolver, bdiag, t, x, lb, ub, (i*bs), bs, constrained_nodal_coords);
                    for (int j=0;j<bs;j++){
                        if(constrained_nodal_coords[j]){
                            ierr = VecSetValue(nbgs->loc_active_indices,(i*bs), 1, INSERT_VALUES); CHKERRQ(ierr);
//                            ierr = VecAssemblyBegin(nbgs->loc_active_indices);CHKERRQ(ierr);
//                            ierr = VecAssemblyEnd(nbgs->loc_active_indices);CHKERRQ(ierr);
                        }
                    }
                }

                ibdiag  += 4;
                bdiag  += 4;

                break;
            case 3:
                v2   = a->a + ii[row+1];
                v3   = a->a + ii[row+2];
                sum1 = b[row];
                sum2 = b[row+1];
                sum3 = b[row+2];
                for (n = 0; n<sz-1; n+=2) {
                    i1    = idx[0];
                    i2    = idx[1];
                    idx  += 2;
                    tmp0  = x[i1];
                    tmp1  = x[i2];
                    sum1 -= v1[0] * tmp0 + v1[1] * tmp1; v1 += 2;
                    sum2 -= v2[0] * tmp0 + v2[1] * tmp1; v2 += 2;
                    sum3 -= v3[0] * tmp0 + v3[1] * tmp1; v3 += 2;
                }

                if (n == sz-1) {
                    tmp0  = x[*idx];
                    sum1 -= v1[0] * tmp0;
                    sum2 -= v2[0] * tmp0;
                    sum3 -= v3[0] * tmp0;
                }
                t[row]   = sum1;
                t[row+1] = sum2;
                t[row+2] = sum3;
                x[row++] = sum1*ibdiag[0] + sum2*ibdiag[3] + sum3*ibdiag[6];
                x[row++] = sum1*ibdiag[1] + sum2*ibdiag[4] + sum3*ibdiag[7];
                x[row++] = sum1*ibdiag[2] + sum2*ibdiag[5] + sum3*ibdiag[8];

                if (nbgs->constrained){
                    nbgs->constrain(&nbgs->blocksolver,  bdiag, t, x, lb, ub, (i*bs), bs, constrained_nodal_coords);
                    for (int j=0;j<bs;j++){
                        if(constrained_nodal_coords[j]){
                            ierr = VecSetValue(nbgs->loc_active_indices,(i*bs +j), 1, INSERT_VALUES); CHKERRQ(ierr);
//                            ierr = VecAssemblyBegin(nbgs->loc_active_indices);CHKERRQ(ierr);
//                            ierr = VecAssemblyEnd(nbgs->loc_active_indices);CHKERRQ(ierr);
                        }
                    }
                }
                ibdiag  += 9;
                bdiag  += 9;
                break;
            case 4:
                v2   = a->a + ii[row+1];
                v3   = a->a + ii[row+2];
                v4   = a->a + ii[row+3];
                sum1 = b[row];
                sum2 = b[row+1];
                sum3 = b[row+2];
                sum4 = b[row+3];
                for (n = 0; n<sz-1; n+=2) {
                    i1    = idx[0];
                    i2    = idx[1];
                    idx  += 2;
                    tmp0  = x[i1];
                    tmp1  = x[i2];
                    sum1 -= v1[0] * tmp0 + v1[1] * tmp1; v1 += 2;
                    sum2 -= v2[0] * tmp0 + v2[1] * tmp1; v2 += 2;
                    sum3 -= v3[0] * tmp0 + v3[1] * tmp1; v3 += 2;
                    sum4 -= v4[0] * tmp0 + v4[1] * tmp1; v4 += 2;
                }

                if (n == sz-1) {
                    tmp0  = x[*idx];
                    sum1 -= v1[0] * tmp0;
                    sum2 -= v2[0] * tmp0;
                    sum3 -= v3[0] * tmp0;
                    sum4 -= v4[0] * tmp0;
                }
                t[row]   = sum1;
                t[row+1] = sum2;
                t[row+2] = sum3;
                t[row+3] = sum4;
                x[row++] = sum1*ibdiag[0] + sum2*ibdiag[4] + sum3*ibdiag[8] + sum4*ibdiag[12];
                x[row++] = sum1*ibdiag[1] + sum2*ibdiag[5] + sum3*ibdiag[9] + sum4*ibdiag[13];
                x[row++] = sum1*ibdiag[2] + sum2*ibdiag[6] + sum3*ibdiag[10] + sum4*ibdiag[14];
                x[row++] = sum1*ibdiag[3] + sum2*ibdiag[7] + sum3*ibdiag[11] + sum4*ibdiag[15];

                if (nbgs->constrained){
                    nbgs->constrain(&nbgs->blocksolver, bdiag, t, x, lb, ub, (i*bs), bs, constrained_nodal_coords);
                    for (int j=0;j<bs;j++){
                        if(constrained_nodal_coords[j]){
                            ierr = VecSetValue(nbgs->loc_active_indices,(i*bs), 1, INSERT_VALUES); CHKERRQ(ierr);
//                            ierr = VecAssemblyBegin(nbgs->loc_active_indices);CHKERRQ(ierr);
//                            ierr = VecAssemblyEnd(nbgs->loc_active_indices);CHKERRQ(ierr);
                        }
                    }
                }
                ibdiag += 16;
                bdiag +=16;
                break;
            case 5:
                v2   = a->a + ii[row+1];
                v3   = a->a + ii[row+2];
                v4   = a->a + ii[row+3];
                v5   = a->a + ii[row+4];
                sum1 = b[row];
                sum2 = b[row+1];
                sum3 = b[row+2];
                sum4 = b[row+3];
                sum5 = b[row+4];
                for (n = 0; n<sz-1; n+=2) {
                    i1    = idx[0];
                    i2    = idx[1];
                    idx  += 2;
                    tmp0  = x[i1];
                    tmp1  = x[i2];
                    sum1 -= v1[0] * tmp0 + v1[1] * tmp1; v1 += 2;
                    sum2 -= v2[0] * tmp0 + v2[1] * tmp1; v2 += 2;
                    sum3 -= v3[0] * tmp0 + v3[1] * tmp1; v3 += 2;
                    sum4 -= v4[0] * tmp0 + v4[1] * tmp1; v4 += 2;
                    sum5 -= v5[0] * tmp0 + v5[1] * tmp1; v5 += 2;
                }

                if (n == sz-1) {
                    tmp0  = x[*idx];
                    sum1 -= v1[0] * tmp0;
                    sum2 -= v2[0] * tmp0;
                    sum3 -= v3[0] * tmp0;
                    sum4 -= v4[0] * tmp0;
                    sum5 -= v5[0] * tmp0;
                }
                t[row]   = sum1;
                t[row+1] = sum2;
                t[row+2] = sum3;
                t[row+3] = sum4;
                t[row+4] = sum5;

                x[row++] = sum1*ibdiag[0] + sum2*ibdiag[5] + sum3*ibdiag[10] + sum4*ibdiag[15] + sum5*ibdiag[20];
                x[row++] = sum1*ibdiag[1] + sum2*ibdiag[6] + sum3*ibdiag[11] + sum4*ibdiag[16] + sum5*ibdiag[21];
                x[row++] = sum1*ibdiag[2] + sum2*ibdiag[7] + sum3*ibdiag[12] + sum4*ibdiag[17] + sum5*ibdiag[22];
                x[row++] = sum1*ibdiag[3] + sum2*ibdiag[8] + sum3*ibdiag[13] + sum4*ibdiag[18] + sum5*ibdiag[23];
                x[row++] = sum1*ibdiag[4] + sum2*ibdiag[9] + sum3*ibdiag[14] + sum4*ibdiag[19] + sum5*ibdiag[24];

                if (nbgs->constrained){
                    nbgs->constrain(&nbgs->blocksolver, bdiag, t, x, lb, ub, (i*bs), bs, constrained_nodal_coords);
                    for (int j=0;j<bs;j++){
                        if(constrained_nodal_coords[j]){
                            ierr = VecSetValue(nbgs->loc_active_indices,(i*bs), 1, INSERT_VALUES); CHKERRQ(ierr);
//                            ierr = VecAssemblyBegin(nbgs->loc_active_indices);CHKERRQ(ierr);
//                            ierr = VecAssemblyEnd(nbgs->loc_active_indices);CHKERRQ(ierr);
                        }
                    }
                }

                ibdiag  += 25;
                bdiag +=25;
                break;

            default:
                SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"Size %D not supported",bs);
        }

        xb   = t;
        ierr = PetscLogFlops(a->nz);CHKERRQ(ierr);
    }

    ierr = VecRestoreArray(xx,&x);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(bb,&b);CHKERRQ(ierr);

    if (nbgs->constrained && nbgs->active_coordinates){
        ierr = VecAssemblyBegin(nbgs->loc_active_indices);CHKERRQ(ierr);
        ierr = VecAssemblyEnd(nbgs->loc_active_indices);CHKERRQ(ierr);
        ierr = VecRestoreLocalVector(nbgs->active_coordinates, nbgs->loc_active_indices); CHKERRQ(ierr);
        //        ierr = VecDestroy(&nbgs->loc_active_indices); CHKERRQ(ierr);
    }

    if (nbgs->constrained) {
        ierr = VecRestoreArray(lblb,&lb);CHKERRQ(ierr);
        ierr = VecRestoreArray(ubub,&ub);CHKERRQ(ierr);
    }

    //    ierr = PetscFree3(a->inode.ibdiag,a->inode.bdiag,a->inode.ssor_work);CHKERRQ(ierr);
    ierr = PetscFree(constrained_nodal_coords); CHKERRQ(ierr);
    PetscFunctionReturn(0);
};

#undef __FUNCT__
#define __FUNCT__ "NBGSStep"
// This does:
//  x = x + \sum_j B_j^1*(b_j- A_j*x_j), under lb <= x <= ub, where j is the number of blocks
PetscErrorCode  NBGSStep(NBGS_CTX * nbgs, Mat A, Vec b, Vec x, Vec _lb, Vec _ub ,PetscInt _blocksize){

    PetscFunctionBegin;
    PetscErrorCode ierr;

    //    Mat loc_A; // holds the local diagonal part of A
    //    Vec res, u, lb_u, ub_u, loc_res, loc_x, loc_u, loc_lb, loc_ub;
    if (nbgs == NULL){
        PetscPrintf(PETSC_COMM_WORLD, "NBGS(...) nbgs struct not set up. Exiting. \n");
        exit(0);
    }
    if (nbgs->initialized == PETSC_FALSE){
        PetscPrintf(PETSC_COMM_WORLD, "NBGS(...) nbgs not initialized. Exiting. \n");
        exit(0);
    }

    if (_lb || _ub){
        ierr = VecSet(nbgs->changed_coordinates, 0);CHKERRQ(ierr);
    }

    { // global communication of error: res = b - A*x
        VecScale(x, -1);
        ierr = MatMultAdd(A, x, b, nbgs->res); CHKERRQ(ierr);
        VecScale(x, -1);
    }

    ierr = VecGetLocalVector(nbgs->res, nbgs->loc_res);CHKERRQ(ierr);
    ierr = VecGetLocalVector(nbgs->u, nbgs->loc_u);CHKERRQ(ierr);

    if (nbgs->constrained){

        ierr = VecWAXPY(nbgs->lb_u, -1, x, _lb); CHKERRQ(ierr);
        ierr = VecGetLocalVector(nbgs->lb_u, nbgs->loc_lb); CHKERRQ(ierr);
        ierr = VecWAXPY(nbgs->ub_u, -1, x, nbgs->ub); CHKERRQ(ierr);
        ierr = VecGetLocalVector(nbgs->ub_u, nbgs->loc_ub); CHKERRQ(ierr);

        local_inode_nbgs(nbgs, nbgs->loc_A, nbgs->loc_res, nbgs->loc_u, nbgs->loc_lb, nbgs->loc_ub, _blocksize, 1);

        ierr = VecSum(nbgs->active_coordinates, &nbgs->num_active_coordinates); CHKERRQ(ierr);
        ierr = VecWAXPY(nbgs->changed_coordinates, -1, nbgs->active_coordinates, nbgs->prev_active_coordinates); CHKERRQ(ierr);
        ierr = VecSum(nbgs->changed_coordinates, &nbgs->num_changed_coordinates); CHKERRQ(ierr);
        ierr = VecCopy(nbgs->active_coordinates, nbgs->prev_active_coordinates); CHKERRQ(ierr);

    }else{
        //        ierr = MatSOR(loc_A,loc_res,1.0,SOR_FORWARD_SWEEP,0,1,1,loc_u);CHKERRQ(ierr);
        local_inode_nbgs(nbgs, nbgs->loc_A, nbgs->loc_res, nbgs->loc_u, NULL, NULL, _blocksize, 1);
    }

    ierr = VecRestoreLocalVector(nbgs->res, nbgs->loc_res); CHKERRQ(ierr);
    ierr = VecRestoreLocalVector(nbgs->u, nbgs->loc_u); CHKERRQ(ierr);

    if (nbgs->constrained){
        ierr = VecRestoreLocalVector(_lb, nbgs->loc_lb); CHKERRQ(ierr);
        ierr = VecRestoreLocalVector(_ub, nbgs->loc_ub); CHKERRQ(ierr);
    }

    if (nbgs->with_ls){
        PetscScalar nom, ctAc;
        {   // determine inactive nodes
            ierr = VecSet(nbgs->mask, 0); CHKERRQ(ierr);
            PetscScalar * active_coord_vals, * mask_vals;
            ierr = VecGetArray(nbgs->active_coordinates, &active_coord_vals); CHKERRQ(ierr);
            ierr = VecGetArray(nbgs->mask, &mask_vals); CHKERRQ(ierr);
            for (int i=0;i<nbgs->n;i++){
                if (active_coord_vals[i] == 0){
                    mask_vals[i] = 1;
                }
            }
            ierr = VecRestoreArray(nbgs->active_coordinates, &active_coord_vals); CHKERRQ(ierr);
            ierr = VecRestoreArray(nbgs->mask, &mask_vals); CHKERRQ(ierr);
        }

        // compute step length: alpha = (res,c) / (c,c)_A
        ierr = VecPointwiseMult(nbgs->c, nbgs->u, nbgs->mask); CHKERRQ(ierr);
        ierr = MatMult(A, nbgs->c, nbgs->Ac); CHKERRQ(ierr);
        ierr = VecDot(nbgs->res, nbgs->c, &nom); CHKERRQ(ierr);
        ierr = VecDot(nbgs->c, nbgs->Ac, &ctAc); CHKERRQ(ierr);

        nbgs->alpha = nom/ctAc;
        // put a max() here for numerical reasons or a warning ...
        if (nbgs->alpha <= 0) {
            PetscPrintf(PETSC_COMM_WORLD, "nbgs: combined correction does not reduce energy! alpha = %2.12e . Setting alpha to 1e-10; \n" , nbgs->alpha);
            nbgs->alpha = 1e-10;
        }
    }else{
        // nothing implemented here yet
        nbgs->alpha = 1;
    }

    if (nbgs->alpha <= 0)
        PetscPrintf(PETSC_COMM_WORLD, "NBGS(...): Are you sure you want omega to be <= zero?  \n");

    ierr = VecAXPY(x,nbgs->alpha,nbgs->u); CHKERRQ(ierr); //x^(k+1) = x^k + omega * B^(-1)*r^k
    ierr  = VecSet(nbgs->active_coordinates,0); CHKERRQ(ierr);

    PetscFunctionReturn(0);
};
#undef __FUNCT__


#undef __FUNCT__
#define __FUNCT__ "NBGSDestroy"
PetscErrorCode  NBGSDestroy(NBGS_CTX * nbgs){
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr = VecDestroy(&nbgs->Ac); CHKERRQ(ierr);
    ierr = VecDestroy(&nbgs->mask); CHKERRQ(ierr);
    ierr = VecDestroy(&nbgs->c); CHKERRQ(ierr);
    ierr = VecDestroy(&nbgs->u); CHKERRQ(ierr);
    ierr = VecDestroy(&nbgs->res); CHKERRQ(ierr);
    ierr = VecDestroy(&nbgs->loc_u); CHKERRQ(ierr);
    ierr = VecDestroy(&nbgs->loc_res); CHKERRQ(ierr);
    ierr = VecDestroy(&nbgs->prev_active_coordinates); CHKERRQ(ierr);
    ierr = VecDestroy(&nbgs->active_coordinates); CHKERRQ(ierr);
    ierr = VecDestroy(&nbgs->changed_coordinates); CHKERRQ(ierr);

    if (nbgs->constrained){
        ierr = VecDestroy(&nbgs->lb_u); CHKERRQ(ierr);
        ierr = VecDestroy(&nbgs->ub_u); CHKERRQ(ierr);
        ierr = VecDestroy(&nbgs->loc_active_indices); CHKERRQ(ierr);
        ierr = VecDestroy(&nbgs->loc_lb); CHKERRQ(ierr);
        ierr = VecDestroy(&nbgs->loc_ub); CHKERRQ(ierr);
        PetscFree(nbgs->blocksolver.diag_constrained);
    }


    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "NBGSCreate"
extern PetscErrorCode NBGSCreate(NBGS_CTX * nbgs, Mat _A, Vec _x, Vec _lb, Vec _ub, PetscInt _blocksize, PetscBool _with_ls, PetscBool _constrained, PetscErrorCode (*constrain)(BLOCKSOLVER * blocksolver, MatScalar * diag, PetscScalar * t, PetscScalar * x,PetscScalar * lb, PetscScalar * ub, int row, int bs, PetscBool * constrained)){

    PetscErrorCode ierr;
    PetscFunctionBegin;

    MPI_Comm_size(PETSC_COMM_WORLD,&nbgs->size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&nbgs->rank);
    nbgs->blocksize = _blocksize;
    nbgs->with_ls = _with_ls;
    nbgs->constrained = _constrained;
    nbgs->lb = _lb;
    nbgs->ub = _ub;
    nbgs->constrain = constrain;

    ierr = MatGetSize(_A, &nbgs->M, &nbgs->N ); CHKERRQ(ierr);
    ierr = MatGetLocalSize(_A, &nbgs->m, &nbgs->n); CHKERRQ(ierr);
    ierr = MatGetDiagonalBlock(_A, &nbgs->loc_A); CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_SELF,&nbgs->loc_res); CHKERRQ(ierr);
    ierr = VecSetSizes(nbgs->loc_res,nbgs->m,PETSC_DECIDE); CHKERRQ(ierr);
    ierr = VecSetFromOptions(nbgs->loc_res); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(nbgs->loc_res); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(nbgs->loc_res); CHKERRQ(ierr);
    ierr = VecDuplicate(nbgs->loc_res, &nbgs->loc_u); CHKERRQ(ierr);

    ierr = VecDuplicate(_x, &nbgs->res); CHKERRQ(ierr);
    ierr = VecDuplicate(_x, &nbgs->u); CHKERRQ(ierr);
    ierr = VecDuplicate(_x, &nbgs->c); CHKERRQ(ierr);
    ierr = VecDuplicate(_x, &nbgs->Ac); CHKERRQ(ierr);
    ierr = VecDuplicate(_x, &nbgs->prev_active_coordinates); CHKERRQ(ierr);
    ierr = VecDuplicate(_x, &nbgs->active_coordinates); CHKERRQ(ierr);
    ierr = VecDuplicate(_x, &nbgs->changed_coordinates); CHKERRQ(ierr);
    ierr = VecDuplicate(_x, &nbgs->mask); CHKERRQ(ierr);

    if (nbgs->constrained){
        ierr = VecDuplicate(_x, &nbgs->lb_u); CHKERRQ(ierr);
        ierr = VecDuplicate(_x, &nbgs->ub_u); CHKERRQ(ierr);
        ierr = VecDuplicate(_x, &nbgs->loc_active_indices); CHKERRQ(ierr);
        ierr = VecDuplicate(nbgs->loc_res, &nbgs->loc_lb); CHKERRQ(ierr);
        ierr = VecDuplicate(nbgs->loc_res, &nbgs->loc_ub); CHKERRQ(ierr);
        PetscMalloc1(_blocksize*_blocksize, &nbgs->blocksolver.diag_constrained);
    }

    ierr = VecSet(nbgs->active_coordinates, 0); CHKERRQ(ierr);
    ierr = VecSet(nbgs->prev_active_coordinates, 0); CHKERRQ(ierr);
    ierr = VecSet(nbgs->changed_coordinates, 0); CHKERRQ(ierr);

    nbgs->initialized = PETSC_TRUE;

    PetscFunctionReturn(0);
}
