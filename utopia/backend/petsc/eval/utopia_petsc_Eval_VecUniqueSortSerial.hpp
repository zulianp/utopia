#ifndef UTOPIA_EVAL_VEC_UNIQUE_SORT_SERIAL_HPP
#define UTOPIA_EVAL_VEC_UNIQUE_SORT_SERIAL_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_petsc_Traits.hpp"
#include "utopia_petsc_Backend.hpp"

#include <cassert>

/*! @file
* Petsc language extensions
*/

namespace utopia 
{
	// This function is not efficient at all
	// it copies values of distributed vector to one processor
	// runs unique search on proc 0
	// removes negative and zero values from the sorted vector
	// final vector values are copied to each processor
    template<class Vector>
    class EvalVecUniqueSortSerial<Vector, PETSC>
    {
        public:
            static void apply(const Wrapper<Vector, 1> &x, Wrapper<Vector, 1> &sorted, const int used_values = -1)  
            { 
                PetscErrorCode ierr = 0;

                PetscInt num_elem = (used_values < 0)? size(x).get(0) : used_values; 

                Vec            V_to_zero;
                VecScatter     ctx; // scatter to all 

                ierr = VecScatterCreateToZero(raw_type(x), &ctx, &V_to_zero);                       assert(ierr == 0);
                ierr = VecScatterBegin(ctx, raw_type(x),V_to_zero, INSERT_VALUES, SCATTER_FORWARD); assert(ierr == 0);
                ierr = VecScatterEnd(ctx, raw_type(x),V_to_zero, INSERT_VALUES, SCATTER_FORWARD);   assert(ierr == 0);

                PetscScalar *v_values; 
                ierr = VecGetArray(V_to_zero, &v_values);                                           assert(ierr == 0);

                PetscInt size; 
                ierr = VecGetSize(V_to_zero, &size);                                                assert(ierr == 0);

                assert(v_values || size == 0);

                if(v_values) {
                    ierr = PetscSortRemoveDupsReal(&size, v_values);                                assert(ierr == 0);
                }

                auto offset = 0; 
                for(auto i=0; i < size; i++)
                {
                    if(v_values[i]>0)
                    {
                        offset = i;
                        break; 
                    }
                }
                
                size -= offset; 
                size = (num_elem < size)? num_elem : size; 

                Vec sorted_zero; 
                MPI_Comm comm = PetscObjectComm((PetscObject) V_to_zero);
                VecCreate(comm, &sorted_zero); 
                VecSetSizes(sorted_zero, size, PETSC_DETERMINE);

                VecType type; 
                VecGetType(V_to_zero, &type); 
                VecSetType(sorted_zero, type);
                VecSet(sorted_zero, 0.0); 
                
                PetscScalar *scalars_sorted; 
                VecGetArray(sorted_zero, &scalars_sorted);              

                for(auto i=0; i < size; i++)
                    scalars_sorted[i] = v_values[offset + i]; 

                VecRestoreArray(sorted_zero, &scalars_sorted); 

                if(!empty(sorted))
                    VecDestroy(&raw_type(sorted)); 

                VecScatter     ctx1; // scatter to all 

                VecScatterCreateToAll(sorted_zero, &ctx1, &raw_type(sorted));
                VecScatterBegin(ctx1, sorted_zero,raw_type(sorted), INSERT_VALUES, SCATTER_FORWARD);
                VecScatterEnd(ctx1, sorted_zero,raw_type(sorted), INSERT_VALUES, SCATTER_FORWARD);

                VecDestroy(&V_to_zero); 
                VecDestroy(&sorted_zero); 
                VecScatterDestroy(&ctx); 
                VecScatterDestroy(&ctx1); 
            }
    };

}

#endif //UTOPIA_EVAL_VEC_UNIQUE_SORT_SERIAL_HPP
