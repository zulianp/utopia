#include "utopia_petsc_Eval_Chop.hpp"
#include "utopia_Tensor.hpp"
#include "utopia_petsc_Matrix.hpp"

namespace utopia {

	template<class Matrix>
    void ChopSmallerThan<Matrix, PETSC>::apply(const Tensor<Matrix, 2> &A, const Scalar &eps)
    {
    
        Scalar    *newVals;
        SizeType       *newCols;
        SizeType       rStart=0, rEnd=0, numRows=0, maxRows=0, r=0, colMax = 0;

        MatGetOwnershipRange(raw_type(A), &rStart, &rEnd);
        for (r = rStart; r < rEnd; ++r) 
        {
            SizeType ncols;

            MatGetRow(raw_type(A), r, &ncols, NULL, NULL);
            colMax = PetscMax(colMax, ncols);
            MatRestoreRow(raw_type(A), r, &ncols, NULL, NULL);
        }
         
        numRows = rEnd - rStart;
        MPIU_Allreduce(&numRows, &maxRows, 1, MPIU_INT, MPI_MAX, PetscObjectComm((PetscObject)raw_type(A)));
        PetscMalloc2(colMax,&newCols,colMax,&newVals);
         
        for (r = rStart; r < rStart+maxRows; ++r) 
        {
            const Scalar *vals;
            const SizeType    *cols;
            SizeType           ncols, newcols, c;

            if (r < rEnd) 
            {
                MatGetRow(raw_type(A), r, &ncols, &cols, &vals);
                for (c = 0; c < ncols; ++c) 
                {
                    newCols[c] = cols[c];
                    newVals[c] = (vals[c]) < eps ? 0.0 : vals[c];
                }
               
                newcols = ncols;
                MatRestoreRow(raw_type(A), r, &ncols, &cols, &vals);
                MatSetValues(raw_type(A), 1, &r, newcols, newCols, newVals, INSERT_VALUES);
            }

            MatAssemblyBegin(raw_type(A), MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(raw_type(A), MAT_FINAL_ASSEMBLY);

        }
           
        PetscFree2(newCols,newVals);
    }

	template<class Matrix>
    void ChopGreaterThan<Matrix, PETSC>::apply(const Tensor<Matrix, 2> &A, const Scalar &eps)
    {
    
        Scalar    *newVals;
        SizeType       *newCols;
        SizeType       rStart=0, rEnd=0, numRows=0, maxRows=0, r=0, colMax = 0;

        MatGetOwnershipRange(raw_type(A), &rStart, &rEnd);
        for (r = rStart; r < rEnd; ++r) 
        {
            SizeType ncols;

            MatGetRow(raw_type(A), r, &ncols, NULL, NULL);
            colMax = PetscMax(colMax, ncols);
            MatRestoreRow(raw_type(A), r, &ncols, NULL, NULL);
        }
         
        numRows = rEnd - rStart;
        MPIU_Allreduce(&numRows, &maxRows, 1, MPIU_INT, MPI_MAX, PetscObjectComm((PetscObject)raw_type(A)));
        PetscMalloc2(colMax,&newCols,colMax,&newVals);
         
        for (r = rStart; r < rStart+maxRows; ++r) 
        {
            const Scalar *vals;
            const SizeType    *cols;
            SizeType           ncols, newcols, c;

            if (r < rEnd) 
            {
                MatGetRow(raw_type(A), r, &ncols, &cols, &vals);
                for (c = 0; c < ncols; ++c) 
                {
                    newCols[c] = cols[c];
                    newVals[c] = (vals[c]) > eps ? 0.0 : vals[c];
                }
               
                newcols = ncols;
                MatRestoreRow(raw_type(A), r, &ncols, &cols, &vals);
                MatSetValues(raw_type(A), 1, &r, newcols, newCols, newVals, INSERT_VALUES);
            }

            MatAssemblyBegin(raw_type(A), MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(raw_type(A), MAT_FINAL_ASSEMBLY);

        }
           
        PetscFree2(newCols,newVals);
    }

    template class ChopSmallerThan<PetscMatrix, PETSC>;
    template class ChopGreaterThan<PetscMatrix, PETSC>;
}
