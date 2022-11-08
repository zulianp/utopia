#include "utopia_libmesh_CreateFE.hpp"

namespace utopia {

    void CreateFE<utopia::libmesh::FunctionSpace, utopia::kokkos::FE<LibMeshScalar_t>>::apply(
        const utopia::libmesh::FunctionSpace &space,
        utopia::kokkos::FE<LibMeshScalar_t> &fe,
        const int degree) {
        assert(false);
        Utopia::Abort();
    }

    void CreateFEOnBoundary<utopia::libmesh::FunctionSpace, utopia::kokkos::FE<LibMeshScalar_t>>::apply(
        const utopia::libmesh::FunctionSpace &space,
        FE &fe,
        int order) {
        assert(false);
        Utopia::Abort();
    }

    void CreateFEOnBoundary<utopia::libmesh::FunctionSpace, utopia::kokkos::FE<LibMeshScalar_t>>::apply(
        const utopia::libmesh::FunctionSpace &space,
        FE &fe,
        const std::string &boundary_name,
        int order) {
        assert(false);
        Utopia::Abort();
    }

    void ConvertField<utopia::Field<utopia::libmesh::FunctionSpace>, LibMeshFEField_t>::apply(const FromField &from,
                                                                                              ToField &to) {
        assert(false);
        Utopia::Abort();
    }

    void ConvertField<LibMeshFEField_t, utopia::Field<utopia::libmesh::FunctionSpace>>::apply(const FromField &from,
                                                                                              ToField &to) {
        assert(false);
        Utopia::Abort();
    }

    void LocalToGlobal<utopia::libmesh::FunctionSpace, LibMeshViewDevice_t, PetscMatrix>::apply(
        const utopia::libmesh::FunctionSpace &space,
        const LibMeshViewDevice_t &element_matrices,
        AssemblyMode mode,
        PetscMatrix &matrix) {
        assert(false);
        Utopia::Abort();
    }

    void LocalToGlobal<utopia::libmesh::FunctionSpace, LibMeshViewDevice_t, PetscVector>::apply(
        const utopia::libmesh::FunctionSpace &space,
        const LibMeshViewDevice_t &element_vectors,
        AssemblyMode mode,
        PetscVector &vector) {
        assert(false);
        Utopia::Abort();
    }

    void LocalToGlobal<utopia::libmesh::FunctionSpace, LibMeshViewDevice_t, PetscVector>::apply(
        const utopia::libmesh::FunctionSpace &space,
        const LibMeshViewDevice_t &element_vectors,
        AssemblyMode mode,
        PetscVector &vector,
        const int n_var) {
        assert(false);
        Utopia::Abort();
    }

    void LocalToGlobal<utopia::libmesh::FunctionSpace, LibMeshViewDevice_t, PetscVector>::side_apply(
        const utopia::libmesh::FunctionSpace &space,
        const LibMeshViewDevice_t &element_vectors,
        AssemblyMode mode,
        PetscVector &vector,
        const std::string &part_name) {
        assert(false);
        Utopia::Abort();
    }

}  // namespace utopia
