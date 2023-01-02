#include "utopia_libmesh_RetroCompatibility.hpp"

#include <libmesh/boundary_info.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_tools.h>

namespace utopia {
    namespace libmesh {

        LibMeshBoundingBox bounding_box(const ::libMesh::MeshBase &mesh) {
#if LIBMESH_VERSION_LESS_THAN(1, 3, 0)
            LibMeshBoundingBox bb = ::libMesh::MeshTools::bounding_box(mesh);
#else
            LibMeshBoundingBox bb = ::libMesh::MeshTools::create_bounding_box(mesh);
#endif
            return bb;
        }

        const ::libMesh::Elem *elem_ptr(const ::libMesh::MeshBase &mesh, const ::libMesh::dof_id_type &idx) {
#if LIBMESH_VERSION_LESS_THAN(1, 6, 0)
            return mesh.elem(idx);
#else
            return mesh.elem_ptr(idx);
#endif
        }

        ::libMesh::Elem *elem_ptr(::libMesh::MeshBase &mesh, const ::libMesh::dof_id_type &idx) {
#if LIBMESH_VERSION_LESS_THAN(1, 6, 0)
            return mesh.elem(idx);
#else
            return mesh.elem_ptr(idx);
#endif
        }

        const ::libMesh::Node *node_ptr(const ::libMesh::MeshBase &mesh, const ::libMesh::dof_id_type &idx) {
#if LIBMESH_VERSION_LESS_THAN(1, 6, 0)
            return &mesh.node(idx);
#else
            return mesh.node_ptr(idx);
#endif
        }

        ::libMesh::Node *node_ptr(::libMesh::MeshBase &mesh, const ::libMesh::dof_id_type &idx) {
#if LIBMESH_VERSION_LESS_THAN(1, 6, 0)
            return &mesh.node(idx);
#else
            return mesh.node_ptr(idx);
#endif
        }

        ::libMesh::boundary_id_type boundary_id(const ::libMesh::BoundaryInfo &b_info,
                                                const ::libMesh::Elem *elem,
                                                const ::libMesh::dof_id_type &idx) {
#if LIBMESH_VERSION_LESS_THAN(1, 6, 0)
            return b_info.boundary_id(elem, idx);
#else
            std::vector< ::libMesh::boundary_id_type> ids;
            b_info.boundary_ids(elem, idx, ids);

            if (ids.empty()) {
                return -1;
            } else {
                return ids[0];
            }
#endif
        }

    }  // namespace libmesh
}  // namespace utopia
