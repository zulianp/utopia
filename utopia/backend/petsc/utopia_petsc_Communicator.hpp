#ifndef UTOPIA_PETSC_COMMUNICATOR_HPP
#define UTOPIA_PETSC_COMMUNICATOR_HPP

#include "utopia_Communicator.hpp"
#include <mpi.h>

namespace utopia {

	class PetscCommunicator final : public MPICommunicator {
	public:
		
		inline Communicator * clone() const override
		{
			return new PetscCommunicator(get());
		}

		inline MPI_Comm get() const override
		{
			return comm_;
		}

		inline void set(MPI_Comm comm)
		{
			comm_ = comm;
		}

		PetscCommunicator(const MPI_Comm comm) : comm_(comm) {}
		PetscCommunicator();

	private:
		MPI_Comm comm_;

	};
}

#endif //UTOPIA_PETSC_COMMUNICATOR_HPP
