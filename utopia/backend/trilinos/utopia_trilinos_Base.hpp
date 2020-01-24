#ifndef UTOPIA_TRILINOS_BASE_HPP
#define UTOPIA_TRILINOS_BASE_HPP

#include "utopia_Base.hpp"
#include <Tpetra_Operator.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <vector>

namespace utopia {

#ifdef UTOPIA_TPETRA_SCALAR
    using TpetraScalar = UTOPIA_TPETRA_SCALAR;
#else
    using TpetraScalar = Tpetra::Operator<>::scalar_type;
#endif

#ifdef UTOPIA_TPETRA_LOCAL_SIZE_TYPE
    using TpetraLocalSizeType = UTOPIA_TPETRA_LOCAL_SIZE_TYPE;
#else
    using TpetraLocalSizeType = Tpetra::Operator<TpetraScalar>::local_ordinal_type;
#endif

#ifdef UTOPIA_TPETRA_SIZE_TYPE
    using TpetraSizeType = UTOPIA_TPETRA_SIZE_TYPE;
#else
    using TpetraSizeType = Tpetra::Operator<TpetraScalar, TpetraLocalSizeType>::global_ordinal_type;
#endif

    //FIXME use Kokkos compatible wrapper
    using TpetraIndexSet      = std::vector<TpetraSizeType>;
    using TpetraIndexArray    = std::vector<TpetraSizeType>;
    using TpetraScalarArray   = std::vector<TpetraScalar>;

    using SerialNode = Kokkos::Compat::KokkosSerialWrapperNode;

#ifdef KOKKOS_ENABLE_CUDA
        using CudaNode   = Kokkos::Compat::KokkosCudaWrapperNode;
        using DefaultKokkosNode = utopia::CudaNode;
#elif defined KOKKOS_ENABLE_ROCM //Kokkos::Compat::KokkosROCmWrapperNode doesn't exist
        using ROCmNode   = Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::ROCm>;
        using DefaultKokkosNode = utopia::ROCmNode;
#elif defined KOKKOS_ENABLE_OPENMP
        using OpenMPNode = Kokkos::Compat::KokkosOpenMPWrapperNode;
        using DefaultKokkosNode = utopia::OpenMPNode;
#else
        using DefaultKokkosNode = utopia::SerialNode;
#endif

}

#endif //UTOPIA_TRILINOS_BASE_HPP
