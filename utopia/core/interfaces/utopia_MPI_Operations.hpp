#ifndef UTOPIA_MPI_OPERATIONS_HPP
#define UTOPIA_MPI_OPERATIONS_HPP

#include "utopia_Operators.hpp"

#ifdef UTOPIA_ENABLE_MPI
#include <mpi.h>
#endif  // UTOPIA_ENABLE_MPI

namespace utopia {

    template <class Op, int Backend>
    class MPIReduceOp {};

#ifdef UTOPIA_ENABLE_MPI

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

#else  // UTOPIA_ENABLE_MPI

#endif  // UTOPIA_ENABLE_MPI

}  // namespace utopia

#endif  // UTOPIA_MPI_OPERATIONS_HPP
