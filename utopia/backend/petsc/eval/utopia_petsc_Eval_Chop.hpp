#ifndef UTOPIA_EVAL_CHOP_HPP
#define UTOPIA_EVAL_CHOP_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_petsc_Traits.hpp"
#include "utopia_petsc_Backend.hpp"

#include <cassert>

/*! @file
* Petsc language extensions
*/

namespace utopia
{


    template<class Matrix>
    class ChopSmallerThan<Matrix, PETSC>
    {
        public:
            static void apply(const Wrapper<Matrix, 2> &A, const double & eps)
            {
            
                PetscScalar    *newVals;
                PetscInt       *newCols;
                PetscInt       rStart, rEnd, numRows, maxRows, r, colMax = 0;

                MatGetOwnershipRange(raw_type(A), &rStart, &rEnd);
                for (r = rStart; r < rEnd; ++r) 
                {
                    PetscInt ncols;

                    MatGetRow(raw_type(A), r, &ncols, NULL, NULL);
                    colMax = PetscMax(colMax, ncols);
                    MatRestoreRow(raw_type(A), r, &ncols, NULL, NULL);
                }
                 
                numRows = rEnd - rStart;
                MPIU_Allreduce(&numRows, &maxRows, 1, MPIU_INT, MPI_MAX, PetscObjectComm((PetscObject)raw_type(A)));
                PetscMalloc2(colMax,&newCols,colMax,&newVals);
                 
                for (r = rStart; r < rStart+maxRows; ++r) 
                {
                    const PetscScalar *vals;
                    const PetscInt    *cols;
                    PetscInt           ncols, newcols, c;

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
    };

    template<class Matrix>
    class ChopBiggerThan<Matrix, PETSC>
    {
        public:
            static void apply(const Wrapper<Matrix, 2> &A, const double & eps)
            {
            
                PetscScalar    *newVals;
                PetscInt       *newCols;
                PetscInt       rStart, rEnd, numRows, maxRows, r, colMax = 0;

                MatGetOwnershipRange(raw_type(A), &rStart, &rEnd);
                for (r = rStart; r < rEnd; ++r) 
                {
                    PetscInt ncols;

                    MatGetRow(raw_type(A), r, &ncols, NULL, NULL);
                    colMax = PetscMax(colMax, ncols);
                    MatRestoreRow(raw_type(A), r, &ncols, NULL, NULL);
                }
                 
                numRows = rEnd - rStart;
                MPIU_Allreduce(&numRows, &maxRows, 1, MPIU_INT, MPI_MAX, PetscObjectComm((PetscObject)raw_type(A)));
                PetscMalloc2(colMax,&newCols,colMax,&newVals);
                 
                for (r = rStart; r < rStart+maxRows; ++r) 
                {
                    const PetscScalar *vals;
                    const PetscInt    *cols;
                    PetscInt           ncols, newcols, c;

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
    };


}

#endif //UTOPIA_EVAL_CHOP_HPP
