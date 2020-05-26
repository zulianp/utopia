#include "utopia_Communicator.hpp"

namespace utopia {

#ifdef WITH_MPI
    bool MPICommunicator::conjunction(const bool &val) const {
        int int_val = val;
        MPI_Allreduce(MPI_IN_PLACE, &int_val, 1, MPI_INT, MPI_SUM, get());

        int size;
        MPI_Comm_size(get(), &size);
        return int_val == size;
    }

    bool MPICommunicator::disjunction(const bool &val) const {
        int int_val = val;
        MPI_Allreduce(MPI_IN_PLACE, &int_val, 1, MPI_INT, MPI_SUM, get());

        int size;
        MPI_Comm_size(get(), &size);
        return int_val > 0;
    }

#endif  // WITH_MPI

}  // namespace utopia