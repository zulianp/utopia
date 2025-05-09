#ifndef UTOPIA_LIBMESH_CREATE_FE_HPP
#define UTOPIA_LIBMESH_CREATE_FE_HPP

#include "utopia_CreateFE.hpp"
#include "utopia_LocalToGlobal.hpp"

#include "utopia_kokkos_FE.hpp"

#include "utopia_libmesh_kokkos_Base.hpp"

namespace utopia {

    template <>
    class CreateFE<utopia::libmesh::FunctionSpace, utopia::kokkos::FE<LibMeshScalar_t>> {
    public:
        using FE = utopia::kokkos::FE<LibMeshScalar_t>;

        static void apply(const utopia::libmesh::FunctionSpace &space,
                          utopia::kokkos::FE<LibMeshScalar_t> &fe,
                          const int order = 0,
                          const int var = 0);
    };

    template <>
    class CreateFEOnBoundary<utopia::libmesh::FunctionSpace, utopia::kokkos::FE<LibMeshScalar_t>> {
    public:
        using FE = utopia::kokkos::FE<LibMeshScalar_t>;

        static void apply(const utopia::libmesh::FunctionSpace &space, FE &fe, int order = 0, int var = 0);
        static void apply(const utopia::libmesh::FunctionSpace &space,
                          FE &fe,
                          const std::string &boundary_name,
                          int order = 0,
                          int var = 0);
    };

    template <>
    class ConvertField<utopia::Field<utopia::libmesh::FunctionSpace>, LibMeshFEField_t> {
    public:
        using FromField = utopia::Field<utopia::libmesh::FunctionSpace>;
        using ToField = LibMeshFEField_t;

        static void apply(const FromField &from, ToField &to);
    };

    template <>
    class ConvertField<LibMeshFEField_t, utopia::Field<utopia::libmesh::FunctionSpace>> {
    public:
        using FromField = LibMeshFEField_t;
        using ToField = utopia::Field<utopia::libmesh::FunctionSpace>;

        static void apply(const FromField &from, ToField &to);
    };

}  // namespace utopia

#endif  // UTOPIA_LIBMESH_CREATE_FE_HPP
