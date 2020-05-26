#ifndef UTOPIA_MPI_OPERATIONS_HPP
#define UTOPIA_MPI_OPERATIONS_HPP

#include "utopia_Operators.hpp"

#ifdef WITH_MPI
#include <mpi.h>
#endif  // WITH_MPI

namespace utopia {

    template <class Op, int Backend>
    class MPIReduceOp {};

#ifdef WITH_MPI

    template <int Backend>
    class MPIReduceOp<Plus, Backend> {
    public:
        inline static MPI_Op op() { return MPI_SUM; }
    };

    template <int Backend>
    class MPIReduceOp<Min, Backend> {
    public:
        inline static MPI_Op op() { return MPI_MIN; }
    };

    template <int Backend>
    class MPIReduceOp<Max, Backend> {
    public:
        inline static MPI_Op op() { return MPI_MAX; }
    };

#else  // WITH_MPI

#endif  // WITH_MPI

}  // namespace utopia

#endif  // UTOPIA_MPI_OPERATIONS_HPP
