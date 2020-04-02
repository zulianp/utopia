#ifndef UTOPIA_PETSC_COMMUNICATOR_HPP
#define UTOPIA_PETSC_COMMUNICATOR_HPP

#include "utopia_Communicator.hpp"
#include <mpi.h>

namespace utopia {

	class PetscCommunicator final : public MPICommunicator {
	public:

		inline PetscCommunicator * clone() const override
		{
			return new PetscCommunicator(get());
		}

		inline MPI_Comm get() const override
		{
			return comm_;
		}

		inline MPI_Comm raw_comm() const override
		{
			return comm_;
		}

		inline void set(MPI_Comm comm)
		{
			comm_ = comm;
		}

		static PetscCommunicator self();

		inline static PetscCommunicator &get_default()
		{
			static PetscCommunicator instance_;
			return instance_;
		}


		explicit PetscCommunicator(const MPI_Comm comm) : comm_(comm) {}
		PetscCommunicator();

		PetscCommunicator(const Communicator &comm) : comm_(comm.raw_comm()) {}

	private:
		MPI_Comm comm_;

	};
}

#endif //UTOPIA_PETSC_COMMUNICATOR_HPP
