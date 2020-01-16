#include "utopia_petsc_Communicator.hpp"

//FIXME
#include "petscvec.h"

namespace utopia {
	PetscCommunicator::PetscCommunicator()
	: comm_(PETSC_COMM_WORLD)
	{

	}
}
