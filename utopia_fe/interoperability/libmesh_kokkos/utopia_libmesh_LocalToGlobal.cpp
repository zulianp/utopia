#include "utopia_libmesh_LocalToGlobal.hpp"

#include "utopia_libmesh_RetroCompatibility.hpp"

// Dofs
#include <libmesh/dof_map.h>

// Mesh
#include <libmesh/boundary_info.h>
#include <libmesh/elem.h>
#include <libmesh/mesh.h>
#include <libmesh/node.h>

namespace utopia {

    void LocalToGlobal<utopia::libmesh::FunctionSpace, LibMeshViewDevice_t, PetscMatrix>::apply(
        const utopia::libmesh::FunctionSpace &space,
        const LibMeshViewDevice_t &d_element_matrices,
        AssemblyMode mode,
        PetscMatrix &matrix,
        int var) {
        // Petsc type defs
        using Scalar = Traits<PetscMatrix>::Scalar;
        using IndexArray = Traits<PetscMatrix>::IndexArray;

        auto &&mesh = space.mesh().raw_type();
        auto &&dof_map = space.raw_type_dof_map();

        const auto n_local_dofs = dof_map.n_local_dofs();
        const long n_cells = mesh.n_active_local_elem();
        auto ebegin = mesh.active_local_elements_begin();
        auto eend = mesh.active_local_elements_end();

        // Scoped lock for parallel assembly
        Write<PetscMatrix> w(matrix, utopia::GLOBAL_ADD);
        if (ebegin == eend) return;

        std::vector<::libMesh::dof_id_type> dof_indices;
        dof_map.dof_indices(*ebegin, dof_indices, var);

        IndexArray idx(dof_indices.size(), 0);
        std::vector<Scalar> values(dof_indices.size() * dof_indices.size(), 0);

        LibMeshViewDevice_t::HostMirror h_element_matrices = Kokkos::create_mirror_view(d_element_matrices);
        Kokkos::deep_copy(h_element_matrices, d_element_matrices);

        int nfun_x_nvar = h_element_matrices.extent(1);
        int n_var = nfun_x_nvar / dof_indices.size();

        assert(var + n_var <= space.n_var());

        SizeType eidx = -1;
        for (auto it = ebegin; it != eend; ++it) {
            auto elem_ptr = *it;
            ++eidx;

            for (int v = 0; v < n_var; ++v) {
                dof_map.dof_indices(*ebegin, dof_indices, var + v);

                int n_shape_functions = dof_indices.size();

                // Convert to PETSc types
                for (int i = 0; i < n_shape_functions; ++i) {
                    idx[i] = dof_indices[i];
                }

                for (int i = 0; i < n_shape_functions; ++i) {
                    for (int j = 0; j < n_shape_functions; ++j) {
                        auto value = h_element_matrices(eidx, i * n_var + v, j * n_var + v);
                        values[i * n_shape_functions + j] = value;
                    }
                }

                matrix.add_matrix(idx, idx, values);
            }
        }
    }

    void LocalToGlobal<utopia::libmesh::FunctionSpace, LibMeshViewDevice_t, PetscVector>::apply(
        const utopia::libmesh::FunctionSpace &space,
        const LibMeshViewDevice_t &d_element_vectors,
        AssemblyMode mode,
        PetscVector &vector,
        int var) {
        apply(space, d_element_vectors, mode, vector, -1, var);
    }

    void LocalToGlobal<utopia::libmesh::FunctionSpace, LibMeshViewDevice_t, PetscVector>::apply(
        const utopia::libmesh::FunctionSpace &space,
        const LibMeshViewDevice_t &d_element_vectors,
        AssemblyMode mode,
        PetscVector &vector,
        const int check_n_var,  // n_var is automatically determined inside!
        int var) {
        // Petsc type defs
        using Scalar = Traits<PetscVector>::Scalar;
        using IndexArray = Traits<PetscVector>::IndexArray;

        auto &&mesh = space.mesh().raw_type();
        auto &&dof_map = space.raw_type_dof_map();

        const auto n_local_dofs = dof_map.n_local_dofs();
        const long n_cells = mesh.n_active_local_elem();
        auto ebegin = mesh.active_local_elements_begin();
        auto eend = mesh.active_local_elements_end();

        // Scoped lock for parallel assembly
        Write<PetscVector> w(vector, utopia::GLOBAL_ADD);
        if (ebegin == eend) return;

        std::vector<::libMesh::dof_id_type> dof_indices;
        dof_map.dof_indices(*ebegin, dof_indices, var);

        IndexArray idx(dof_indices.size(), 0);
        std::vector<Scalar> values(dof_indices.size(), 0);

        LibMeshViewDevice_t::HostMirror h_element_vectors = Kokkos::create_mirror_view(d_element_vectors);
        Kokkos::deep_copy(h_element_vectors, d_element_vectors);

        int nfun_x_nvar = h_element_vectors.extent(1);
        int n_var = nfun_x_nvar / dof_indices.size();

        assert(check_n_var == -1 || n_var == check_n_var);

        assert(var + n_var <= space.n_var());

        SizeType eidx = -1;
        for (auto it = ebegin; it != eend; ++it) {
            auto elem_ptr = *it;
            ++eidx;

            for (int v = 0; v < n_var; ++v) {
                dof_map.dof_indices(elem_ptr, dof_indices, var + v);

                int n_shape_functions = dof_indices.size();

                // Convert to PETSc types
                for (int i = 0; i < n_shape_functions; ++i) {
                    idx[i] = dof_indices[i];
                }

                for (int i = 0; i < n_shape_functions; ++i) {
                    auto value = h_element_vectors(eidx, i * n_var + v);
                    values[i] = value;
                }

                vector.add_vector(idx, values);
            }
        }
    }

    void LocalToGlobal<utopia::libmesh::FunctionSpace, LibMeshViewDevice_t, PetscVector>::side_apply(
        const utopia::libmesh::FunctionSpace &space,
        const LibMeshViewDevice_t &d_element_vectors,
        AssemblyMode mode,
        PetscVector &vector,
        const std::string &boundary_name,
        int var) {
        // Petsc type defs
        using Scalar = Traits<PetscVector>::Scalar;
        using IndexArray = Traits<PetscVector>::IndexArray;

        auto &&mesh = space.mesh().raw_type();
        auto &&dof_map = space.raw_type_dof_map();

        const auto n_local_dofs = dof_map.n_local_dofs();
        const long n_cells = mesh.n_active_local_elem();
        auto ebegin = mesh.active_local_elements_begin();
        auto eend = mesh.active_local_elements_end();

        auto &&binfo = mesh.get_boundary_info();

        int boundary_id_select = -1;
        if (!boundary_name.empty()) {
            boundary_id_select = binfo.get_id_by_name(boundary_name);
        }

        // Scoped lock for parallel assembly
        Write<PetscVector> w(vector, utopia::GLOBAL_ADD);
        if (ebegin == eend) return;

        libMesh::Elem *boundary_element = nullptr;
        int boundary_side = 0;

        {
            // Find one boundary element!
            for (const auto &elem_ptr : mesh.active_local_element_ptr_range()) {
                const int n_sides = elem_ptr->n_sides();

                for (int i = 0; i < n_sides; ++i) {
                    if ((elem_ptr->neighbor_ptr(i) != libmesh_nullptr)) {
                        continue;
                    }

                    bool selected = boundary_id_select == -1 ||
                                    boundary_id_select == utopia::libmesh::boundary_id(binfo, elem_ptr, i);

                    if (!selected) continue;

                    if (!boundary_element) {
                        // We have found a boundary element!
                        boundary_element = elem_ptr;
                        boundary_side = i;
                        break;
                    }
                }

                if (boundary_element) break;
            }
        }

        if (!boundary_element) return;

        // std::unique_ptr<libMesh::Elem> side_ptr = boundary_element->build_side_ptr(boundary_side);

        // std::vector<::libMesh::dof_id_type> dof_indices;
        // dof_map.dof_indices(side_ptr, dof_indices, var);

        // IndexArray idx(dof_indices.size(), 0);
        // std::vector<Scalar> values(dof_indices.size(), 0);

        // LibMeshViewDevice_t::HostMirror h_element_vectors = Kokkos::create_mirror_view(d_element_vectors);
        // Kokkos::deep_copy(h_element_vectors, d_element_vectors);

        // int nfun_x_nvar = h_element_vectors.extent(1);
        // int n_var = nfun_x_nvar / dof_indices.size();

        // assert(var + n_var <= space.n_var());

        // SizeType eidx = 0;
        // for (auto it = ebegin; it != eend; ++it) {
        //     auto elem_ptr = *it;

        //     for (int i = 0; i < n_sides; ++i) {
        //         if ((elem_ptr->neighbor_ptr(i) != libmesh_nullptr)) {
        //             continue;
        //         }

        //         bool selected =
        //             boundary_id_select == -1 || boundary_id_select == utopia::libmesh::boundary_id(binfo, elem_ptr,
        //             i);

        //         if (!selected) continue;

        //         std::unique_ptr<libMesh::Elem> side_ptr = elem_ptr->build_side_ptr(i);

        //         for (int v = 0; v < n_var; ++v) {
        //             dof_map.dof_indices(side_ptr.get(), dof_indices, var + v);

        //             int n_shape_functions = dof_indices.size();

        //             // Convert to PETSc types
        //             for (int i = 0; i < n_shape_functions; ++i) {
        //                 idx[i] = dof_indices[i];
        //             }

        //             for (int i = 0; i < n_shape_functions; ++i) {
        //                 auto value = h_element_vectors(eidx, i * n_var + v);
        //                 values[i] = value;
        //             }

        //             vector.add_vector(idx, values);
        //             eidx++;
        //         }
        //     }
        // }
    }

}  // namespace utopia
