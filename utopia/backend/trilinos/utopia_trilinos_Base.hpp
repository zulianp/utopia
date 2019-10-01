#ifndef UTOPIA_TRILINOS_BASE_HPP
#define UTOPIA_TRILINOS_BASE_HPP

#include <Tpetra_Operator.hpp>
#include <vector>

#define UTOPIA_LAMBDA KOKKOS_LAMBDA

namespace utopia {

    using TpetraScalar        = Tpetra::Operator<>::scalar_type;
    using TpetraLocalSizeType = Tpetra::Operator<TpetraScalar>::local_ordinal_type;
    using TpetraSizeType      = Tpetra::Operator<TpetraScalar, TpetraLocalSizeType>::global_ordinal_type;

    using TpetraIndexSet      = std::vector<TpetraSizeType>;
    using TpetraIndexArray    = std::vector<TpetraSizeType>;
    using TpetraScalarArray   = std::vector<TpetraScalar>;

}

#endif //UTOPIA_TRILINOS_BASE_HPP
