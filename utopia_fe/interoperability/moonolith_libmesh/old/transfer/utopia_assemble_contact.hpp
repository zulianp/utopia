#ifndef UTOPIA_ASSEMBLE_CONTACT_HPP
#define UTOPIA_ASSEMBLE_CONTACT_HPP

#include <memory>
#include <utility>
#include <vector>

#include "libmesh/libmesh_common.h"
#include "utopia.hpp"
#include "utopia_fe_base.hpp"

// forward decl
namespace moonolith {
    class Communicator;
}

namespace libMesh {
    class MeshBase;
    class DofMap;
}  // namespace libMesh

namespace utopia {

    bool assemble_contact(moonolith::Communicator &comm,
                          const std::shared_ptr<libMesh::MeshBase> &mesh,
                          const std::shared_ptr<libMesh::DofMap> &dof_map,
                          const unsigned int var_num,
                          USparseMatrix &B,
                          USparseMatrix &orthogonal_trafos,
                          UVector &gap,
                          USparseMatrix &normals,
                          UVector &is_contact_node,
                          const libMesh::Real search_radius,
                          const std::vector<std::pair<int, int> > &tags,
                          const bool use_biorth = true,
                          const bool use_volume_differential = false);

    bool assemble_contact(moonolith::Communicator &comm,
                          const std::shared_ptr<libMesh::MeshBase> &mesh,
                          const std::shared_ptr<libMesh::DofMap> &dof_map,
                          const unsigned int var_num,
                          USparseMatrix &B,
                          USparseMatrix &orthogonal_trafos,
                          UVector &gap,
                          USparseMatrix &normals,
                          UVector &is_contact_node,
                          const libMesh::Real search_radius,
                          const int tag_1,
                          const int tag_2,
                          const bool use_biorth = true,
                          const bool use_volume_differential = false);

    bool assemble_contact(moonolith::Communicator &comm,
                          const std::shared_ptr<libMesh::MeshBase> &mesh,
                          const std::shared_ptr<libMesh::DofMap> &dof_map,
                          const unsigned int var_num,
                          USparseMatrix &B,
                          USparseMatrix &orthogonal_trafos,
                          UVector &gap,
                          UVector &normals,
                          UVector &is_contact_node,
                          const libMesh::Real search_radius,
                          const std::vector<std::pair<int, int> > &tags,
                          const bool use_biorth = true,
                          const bool use_volume_differential = false);

    void convert_normal_matrix_to_vector(const USparseMatrix &mat, UVector &vec);

}  // namespace utopia

#endif  // UTOPIA_ASSEMBLE_CONTACT_HPP
