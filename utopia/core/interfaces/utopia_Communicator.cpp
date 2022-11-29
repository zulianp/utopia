#include "utopia_Communicator.hpp"

namespace utopia {

#ifdef UTOPIA_WITH_MPI
    bool MPIBaseCommunicator::conjunction(const bool &val) const {
        int int_val = val;
        MPI_Allreduce(MPI_IN_PLACE, &int_val, 1, MPI_INT, MPI_SUM, get());

        int size;
        MPI_Comm_size(get(), &size);
        return int_val == size;
    }

    bool MPIBaseCommunicator::disjunction(const bool &val) const {
        int int_val = val;
        MPI_Allreduce(MPI_IN_PLACE, &int_val, 1, MPI_INT, MPI_SUM, get());

        int size;
        MPI_Comm_size(get(), &size);
        return int_val > 0;
    }

    MPICommunicator MPICommunicator::split(const int color) const {
        MPI_Comm row_comm;
        MPI_Comm_split(get(), color, rank(), &row_comm);
        MPICommunicator ret(row_comm, true);
        return ret;
    }

#endif  // UTOPIA_WITH_MPI

}  // namespace utopia