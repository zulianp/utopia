#ifndef UTOPIA_MPI_HPP
#define UTOPIA_MPI_HPP 

#include "utopia_Base.hpp"

namespace utopia {
	SizeType mpi_world_size();
	SizeType mpi_world_rank();
	void mpi_world_barrier();
}

#endif //UTOPIA_MPI_HPP
