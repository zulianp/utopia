#ifndef UTOPIA_LIBMESH_LOCAL_TO_GLOBAL_HPP
#define UTOPIA_LIBMESH_LOCAL_TO_GLOBAL_HPP

#include "utopia_LocalToGlobal.hpp"

#include "utopia_kokkos_FE.hpp"

#include "utopia_libmesh_kokkos_Base.hpp"

namespace utopia {

    template <>
    class LocalToGlobal<utopia::libmesh::FunctionSpace, LibMeshViewDevice_t, PetscMatrix> {
    public:
        static void apply(const utopia::libmesh::FunctionSpace &space,
                          const LibMeshViewDevice_t &element_matrices,
                          AssemblyMode mode,
                          PetscMatrix &matrix,
                          int var = 0);
    };

    template <>
    class LocalToGlobal<utopia::libmesh::FunctionSpace, LibMeshViewDevice_t, PetscVector> {
    public:
        static void apply(const utopia::libmesh::FunctionSpace &space,
                          const LibMeshViewDevice_t &element_vectors,
                          AssemblyMode mode,
                          PetscVector &vector,
                          int var = 0);

        static void apply(const utopia::libmesh::FunctionSpace &space,
                          const LibMeshViewDevice_t &element_vectors,
                          AssemblyMode mode,
                          PetscVector &vector,
                          const int n_var,
                          int var = 0);

        static void side_apply(const utopia::libmesh::FunctionSpace &space,
                               const LibMeshViewDevice_t &element_vectors,
                               AssemblyMode mode,
                               PetscVector &vector,
                               const std::string &part_name,
                               int var = 0);
    };

}  // namespace utopia

#endif  // UTOPIA_LIBMESH_LOCAL_TO_GLOBAL_HPP
