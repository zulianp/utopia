#ifndef UTOPIA_MAT_LINEAR_SOLVER_HPP
#define UTOPIA_MAT_LINEAR_SOLVER_HPP

#include "utopia_LinearSolver.hpp"
#include "utopia_Core.hpp"

namespace utopia 
{

    // this class needs some refactoring.... and more thinking... 
    template<typename Matrix, typename DenseMatrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class MatLinearSolver {};

    template<typename Matrix, typename DenseMatrix, typename Vector>
    class MatLinearSolver<Matrix, DenseMatrix, Vector, PETSC> 
    {
    public:
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;

        static_assert(Traits<Matrix>::Backend == utopia::PETSC, "MatLinearSolver only works with petsc types");
        static_assert(!utopia::is_sparse<DenseMatrix>::value, "MatLinearSolver does not support sparse matrices for X and RHS."); 

        MatLinearSolver(const std::shared_ptr <LinearSolver> &linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector> >())
        :linear_solver_(linear_solver)
        {

        }

        virtual bool solve(const Matrix &A, const DenseMatrix & B, DenseMatrix & X)
        {
            Vector       b,x;

            PetscInt       m,N,i;
            PetscScalar    *bb,*xx;
            PetscBool      flg;

            PetscObjectTypeCompareAny((PetscObject)raw_type(B),&flg,MATSEQDENSE,MATMPIDENSE,NULL);
            if (!flg) SETERRQ(PetscObjectComm((PetscObject)raw_type(B)),PETSC_ERR_ARG_WRONG,"Matrix B must be MATDENSE matrix");
            PetscObjectTypeCompareAny((PetscObject)raw_type(X),&flg,MATSEQDENSE,MATMPIDENSE,NULL);
            if (!flg) SETERRQ(PetscObjectComm((PetscObject)raw_type(B)),PETSC_ERR_ARG_WRONG,"Matrix X must be MATDENSE matrix");

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

            return false; 
        }




    protected:
        std::shared_ptr<LinearSolver> linear_solver_;
    };

}

#endif //UTOPIA_MAT_LINEAR_SOLVER_HPP
