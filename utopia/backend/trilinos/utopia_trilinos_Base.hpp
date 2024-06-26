#ifndef UTOPIA_TRILINOS_BASE_HPP
#define UTOPIA_TRILINOS_BASE_HPP

#include <Trilinos_version.h>
#if (TRILINOS_MAJOR_MINOR_VERSION >= 140000)
#include <Tpetra_KokkosCompat_DefaultNode.hpp>
#else
#if (TRILINOS_MAJOR_MINOR_VERSION >= 130500)
#include <KokkosCompat_DefaultNode.hpp>
#else
#include <Kokkos_DefaultNode.hpp>
#endif
#endif

#include <Tpetra_Operator.hpp>
#include <vector>
#include "utopia_Base.hpp"

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

    // FIXME use Kokkos compatible wrapper
    using TpetraIndexSet = std::vector<TpetraSizeType>;
    using TpetraIndexArray = std::vector<TpetraSizeType>;
    using TpetraScalarArray = std::vector<TpetraScalar>;

#if (TRILINOS_MAJOR_MINOR_VERSION >= 140000)
    using SerialNode = Tpetra::KokkosCompat::KokkosSerialWrapperNode;
#else
    using SerialNode = Kokkos::Compat::KokkosSerialWrapperNode;
#endif

#ifdef KOKKOS_ENABLE_CUDA
#if (TRILINOS_MAJOR_MINOR_VERSION >= 140000)
    using CudaNode = Tpetra::KokkosCompat::KokkosCudaWrapperNode;
#else
    using CudaNode = Kokkos::Compat::KokkosCudaWrapperNode;
#endif
    using DefaultKokkosNode = utopia::CudaNode;
#elif defined KOKKOS_ENABLE_ROCM  // Kokkos::Compat::KokkosROCmWrapperNode doesn't exist
#if (TRILINOS_MAJOR_MINOR_VERSION >= 140000)
    using ROCmNode = Tpetra::KokkosCompat::KokkosDeviceWrapperNode<Kokkos::ROCm>;
#else
    using ROCmNode = Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::ROCm>;
#endif
    using DefaultKokkosNode = utopia::ROCmNode;
#elif defined KOKKOS_ENABLE_OPENMP
#if (TRILINOS_MAJOR_MINOR_VERSION >= 140000)
    using OpenMPNode = Tpetra::KokkosCompat::KokkosOpenMPWrapperNode;
#else
    using OpenMPNode = Kokkos::Compat::KokkosOpenMPWrapperNode;
#endif
    using DefaultKokkosNode = utopia::OpenMPNode;
#else
    using DefaultKokkosNode = utopia::SerialNode;
#endif

}  // namespace utopia

#endif  // UTOPIA_TRILINOS_BASE_HPP
