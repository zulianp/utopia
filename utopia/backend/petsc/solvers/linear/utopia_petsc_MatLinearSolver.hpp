#ifndef UTOPIA_MAT_LINEAR_SOLVER_HPP
#define UTOPIA_MAT_LINEAR_SOLVER_HPP

#include "utopia_LinearSolver.hpp"
#include "utopia_Core.hpp"

namespace utopia 
{

    template<typename Matrix, typename DenseMatrix, typename Vector>
    class MatLinearSolver<Matrix, DenseMatrix, Vector, PETSC> 
    {
    public:
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;
        static_assert(Traits<Matrix>::Backend == utopia::PETSC, "MatLinearSolver only works with petsc types");

        MatLinearSolver(const std::shared_ptr <LinearSolver> &linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector> >())
        :linear_solver_(linear_solver)
        {

        }

        virtual bool solve(const Matrix &A, const DenseMatrix & B, DenseMatrix & X)
        {
            Vector       b,x;

            PetscInt       m,N,i;
            PetscScalar    *bb,*xx;
            PetscBool      flgX, flgB;

            PetscObjectTypeCompareAny((PetscObject)raw_type(B),&flgB,MATSEQDENSE,MATMPIDENSE,NULL);
            PetscObjectTypeCompareAny((PetscObject)raw_type(X),&flgX,MATSEQDENSE,MATMPIDENSE,NULL);

            if(flgB && flgX)
            {
                MatDenseGetArray(raw_type(B),&bb);
                MatDenseGetArray(raw_type(X),&xx);
                MatGetLocalSize(raw_type(B),&m,NULL);  /* number local rows */
                MatGetSize(raw_type(B),NULL, &N);       /* total columns in dense matrix */
                MatCreateVecs(raw_type(A), &raw_type(x), &raw_type(b));
                for (i=0; i<N; i++) 
                {
                    VecPlaceArray(raw_type(b), bb + i*m);
                    VecPlaceArray(raw_type(x), xx + i*m);
                    linear_solver_->solve(A, b, x); 
                    VecResetArray(raw_type(x));
                    VecResetArray(raw_type(b));
                }
                VecDestroy(&raw_type(b));
                VecDestroy(&raw_type(x));
                MatDenseRestoreArray(raw_type(B), &bb);
                MatDenseRestoreArray(raw_type(X), &xx);
            }
            else
            {
                MatType typeX; 
                MatGetType(raw_type(X), &typeX); 

                Mat B_temp, X_temp; 
                MatConvert(raw_type(B), MATDENSE, MAT_INITIAL_MATRIX, &B_temp); 
                MatConvert(raw_type(X), MATDENSE, MAT_INITIAL_MATRIX, &X_temp); 

                MatDenseGetArray(B_temp,&bb);
                MatDenseGetArray(X_temp,&xx);
                MatGetLocalSize(B_temp,&m,NULL);  /* number local rows */
                MatGetSize(B_temp,NULL, &N);       /* total columns in dense matrix */
                MatCreateVecs(raw_type(A), &raw_type(x), &raw_type(b));

                for (i=0; i<N; i++) 
                {
                    VecPlaceArray(raw_type(b), bb + i*m);
                    VecPlaceArray(raw_type(x), xx + i*m);
                    linear_solver_->solve(A, b, x); 
                    VecResetArray(raw_type(x));
                    VecResetArray(raw_type(b));
                }
                VecDestroy(&raw_type(b));
                VecDestroy(&raw_type(x));
                MatDenseRestoreArray(B_temp, &bb);
                MatDenseRestoreArray(X_temp, &xx);


                MatConvert(X_temp, typeX, MAT_REUSE_MATRIX , &raw_type(X)); 
                
                MatDestroy(&B_temp);
                MatDestroy(&X_temp);         
            }

            return false; 
        }


        // very unefficient routine, precision of the inverse depends on LS
        virtual void get_inverse(const Matrix &A, DenseMatrix & A_inv)
        {
            if(size(A).get(0) != size(A).get(1))
                utopia_error("get_inverse:: works just with square matrices \n"); 

            DenseMatrix RHS = local_identity(local_size(A).get(0), local_size(A).get(1)); 
            A_inv = local_values(local_size(A).get(0), local_size(A).get(1), 0.0);
            this->solve(A, RHS, A_inv); 
        }


    protected:
        std::shared_ptr<LinearSolver> linear_solver_;
    };

}

#endif //UTOPIA_MAT_LINEAR_SOLVER_HPP
