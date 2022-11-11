#ifndef UTOPIA_LIBMESH_RETROCOMPATIBILITY_HPP
#define UTOPIA_LIBMESH_RETROCOMPATIBILITY_HPP

#include "utopia_Base.hpp"

#include "libmesh/bounding_box.h"
#include "libmesh/libmesh_version.h"

// Forward declarations
namespace libMesh {
    class MeshBase;
    class BoundaryInfo;
    class Elem;
    class Node;
}  // namespace libMesh

namespace utopia {
    namespace libmesh {

#if LIBMESH_VERSION_LESS_THAN(1, 6, 0)
        using LibMeshBoundingBox = ::libMesh::MeshTools::BoundingBox;
#else
        using LibMeshBoundingBox = ::libMesh::BoundingBox;
#endif

        LibMeshBoundingBox bounding_box(const ::libMesh::MeshBase &mesh);

        const ::libMesh::Elem *elem_ptr(const ::libMesh::MeshBase &mesh, const ::libMesh::dof_id_type &idx);

        ::libMesh::Elem *elem_ptr(::libMesh::MeshBase &mesh, const ::libMesh::dof_id_type &idx);

        const ::libMesh::Node *node_ptr(const ::libMesh::MeshBase &mesh, const ::libMesh::dof_id_type &idx);

        ::libMesh::Node *node_ptr(::libMesh::MeshBase &mesh, const ::libMesh::dof_id_type &idx);

        UTOPIA_DEPRECATED_MSG("This method should be replaced. multiple ids have to be handled.")
        ::libMesh::boundary_id_type boundary_id(const ::libMesh::BoundaryInfo &b_info,
                                                const ::libMesh::Elem *elem,
                                                const ::libMesh::dof_id_type &idx);

    }  // namespace libmesh
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_RETROCOMPATIBILITY_HPP
