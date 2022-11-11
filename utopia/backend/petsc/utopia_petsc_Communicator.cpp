#include "utopia_petsc_Communicator.hpp"

// FIXME
#include "petscvec.h"

namespace utopia {

    PetscCommunicator::PetscCommunicator() : MPICommunicator(PETSC_COMM_WORLD) {}

    PetscCommunicator PetscCommunicator::self() { return PetscCommunicator(PETSC_COMM_SELF); }

    PetscCommunicator PetscCommunicator::world() { return PetscCommunicator(PETSC_COMM_WORLD); }

}  // namespace utopia
