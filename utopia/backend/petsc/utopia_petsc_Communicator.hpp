#ifndef UTOPIA_PETSC_COMMUNICATOR_HPP
#define UTOPIA_PETSC_COMMUNICATOR_HPP

#include "utopia_Communicator.hpp"
#include <mpi.h>

namespace utopia {

	class PetscCommunicator final : public Communicator {
	public:

		inline int rank() const override
		{
			int ret;
			MPI_Comm_rank(get(), &ret);
			return ret;
		}

		inline int size() const override
		{
			int ret;
			MPI_Comm_size(get(), &ret);
			return ret;
		}

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
