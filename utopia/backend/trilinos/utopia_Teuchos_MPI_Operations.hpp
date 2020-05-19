#ifndef UTOPIA_TEUCHOS_MPI_OPERATIONS_HPP
#define UTOPIA_TEUCHOS_MPI_OPERATIONS_HPP

#include "utopia_MPI_Operations.hpp"
#include "utopia_trilinos_Base.hpp"

#include <Teuchos_EReductionType.hpp>

namespace utopia {

    template <>
    class MPIReduceOp<Plus, TRILINOS> {
    public:
        inline static Teuchos::EReductionType op() { return Teuchos::REDUCE_SUM; }
    };

    template <>
    class MPIReduceOp<Min, TRILINOS> {
    public:
        inline static Teuchos::EReductionType op() { return Teuchos::REDUCE_MIN; }
    };

    template <>
    class MPIReduceOp<Max, TRILINOS> {
    public:
        inline static Teuchos::EReductionType op() { return Teuchos::REDUCE_MAX; }
    };
}  // namespace utopia

#endif  // UTOPIA_TEUCHOS_MPI_OPERATIONS_HPP
