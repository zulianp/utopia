#include "utopia_petsc_Vector.hpp"

namespace utopia {

	void PETScVector::repurpose(
			MPI_Comm comm,
			VecType type,
			PetscInt n_local,
			PetscInt n_global)
	{
		if(vec_ == nullptr) {
			VecCreate(comm, &vec_);
		} else {
			if(comm != PetscObjectComm((PetscObject)vec_)) {
				check_error( VecDestroy(&vec_) );
				check_error( VecCreate(comm, &vec_) );
			} else {
				PetscInt old_n_global;
				VecGetSize(vec_, &old_n_global);

				if(old_n_global == n_global) {
					PetscInt old_n_local;
					VecGetLocalSize(vec_, &old_n_local);
					assert(old_n_local == n_local && "We do not handle the local consistency. Explicitly set sizes in the initialization.");
					return;
				}
			}
		}
		
		check_error( VecSetType(vec_, type) );
		check_error( VecSetSizes(vec_, n_local, n_global) );

		ghost_values_.clear();
		
		initialized_ = true;
	}

	void PETScVector::init(MPI_Comm comm, VecType type, PetscInt n_local, PetscInt n_global)
	{
        check_error( VecCreate(comm, &vec_) );
        check_error( VecSetType(vec_, type) );
        check_error( VecSetSizes(vec_, n_local, n_global) );
        
        initialized_ = true;
	}

	void PETScVector::ghosted(
		MPI_Comm comm,
	    PetscInt local_size,
	    PetscInt global_size,
	    const std::vector<PetscInt> &index)
	{
	   	destroy();

	    check_error(
	        VecCreateGhost(
	            comm, 
	            local_size,
	            global_size,
	            static_cast<PetscInt>(index.size()),
	            &index[0],
	            &vec_)
	    );

	    init_ghost_index(index);
	    VecZeroEntries(vec_);

	    initialized_ = true;
	}
}
