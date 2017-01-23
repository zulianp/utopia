#include "utopia_MPI.hpp"


#ifdef WITH_MPI

#include <mpi.h>

namespace utopia 
{
    /** \addtogroup MPI
     * @brief MPI based routines. 
     * @ingroup parallel_expressions
     */

	/**
	 * @ingroup 	MPI
	 * @brief      	Getter for MPI Comm size. 
	 *
	 * @return     MPI_COMM_WORLD size. 
	 */
	SizeType mpi_world_size()
	{
		int size;
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		return size;
	}


	/**
	 * @ingroup 	MPI
	 * @brief      	Getter for MPI world rank.
	 *
	 * @return     The rank on given processor. 
	 */
	SizeType mpi_world_rank()
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		return rank;
	}

	
	/**
	 * @ingroup 	MPI
	 * @brief      	Blocks until all processes in the MPI_COMM_WORLD communicator have reached this routine. 
	 *
	 */
	void mpi_world_barrier()
	{
		MPI_Barrier(MPI_COMM_WORLD);
	}
}

#else
namespace utopia 
{

	SizeType mpi_world_size()
	{
		return 1;
	}

	SizeType mpi_world_rank()
	{
		return 0;
	}

	void mpi_world_barrier()
	{
	
	}
}

#endif //WITH_MPI



