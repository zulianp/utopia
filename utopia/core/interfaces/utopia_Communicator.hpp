#ifndef UTOPIA_COMMUNICATOR_HPP
#define UTOPIA_COMMUNICATOR_HPP

#include "utopia_Base.hpp"

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace utopia {

	class Communicator {
	public:
		virtual ~Communicator() {}
		virtual int rank() const = 0;
		virtual int size() const = 0;
		virtual Communicator * clone() const = 0;

// #ifdef WITH_MPI
// 		virtual MPI_Comm get() const = 0;
// #endif //WITH_MPI

	};

	/**
	 * @brief usefull for non-mpi backends as a default object to return/provide
	 */
	class SelfCommunicator final : Communicator {
	public:
		int rank() const noexcept override { return 0; }
		int size() const noexcept override { return 1; }

		SelfCommunicator * clone() const noexcept override 
		{
			return new SelfCommunicator();
		}

#ifdef WITH_MPI
		inline MPI_Comm get() const //override 
		{
			return MPI_COMM_SELF;
		}
#endif //WITH_MPI
	};


}

#endif //UTOPIA_COMMUNICATOR_HPP
