#include "utopia_petsc_Communicator.hpp"

//FIXME
#include "petscvec.h"

namespace utopia {

	PetscCommunicator::PetscCommunicator()
	: wrapper_( make_not_owned(PETSC_COMM_WORLD) )
	{}

    PetscCommunicator PetscCommunicator::self()
    {
        return PetscCommunicator(PETSC_COMM_SELF);
    }

    PetscCommunicator PetscCommunicator::world()
    {
        return PetscCommunicator(PETSC_COMM_WORLD);
    }

    PetscCommunicator PetscCommunicator::split(const int color) const
    {
        MPI_Comm row_comm;
        MPI_Comm_split(get(), color, rank(), &row_comm);
        PetscCommunicator ret;
        ret.own(row_comm);
        return ret;
    }

}
