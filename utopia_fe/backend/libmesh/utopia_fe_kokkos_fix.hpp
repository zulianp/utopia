#ifndef UTOPIA_FE_KOKKOS_FIX_HPP
#define UTOPIA_FE_KOKKOS_FIX_HPP

#include "utopia_fe_config.hpp"

#ifdef WITH_TRILINOS
// bug in Kokkos-Kernels, this has to be included first (https://github.com/kokkos/kokkos-kernels/issues/309)
#include <Eigen/Core>

#endif //WITH_TRILINOS

#endif //UTOPIA_FE_KOKKOS_FIX_HPP