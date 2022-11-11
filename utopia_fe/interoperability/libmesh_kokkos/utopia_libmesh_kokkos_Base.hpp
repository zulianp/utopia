#ifndef UTOPIA_LIBMESH_KOKKOS_BASE_HPP
#define UTOPIA_LIBMESH_KOKKOS_BASE_HPP

#include "utopia_fe_Core.hpp"

#include "utopia_libmesh_FunctionSpace_new.hpp"

#include "utopia_kokkos_FEBase.hpp"
#include "utopia_kokkos_Field.hpp"

namespace utopia {
    using LibMeshScalar_t = utopia::Traits<utopia::libmesh::FunctionSpace>::Scalar;
    using LibMeshViewDevice_t = utopia::kokkos::DefaultView<LibMeshScalar_t>;
    using LibMeshIntViewDevice_t = utopia::kokkos::DefaultView<int>;
    using LibMeshFEField_t = utopia::kokkos::Field<utopia::kokkos::FE<LibMeshScalar_t>>;
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_KOKKOS_BASE_HPP
