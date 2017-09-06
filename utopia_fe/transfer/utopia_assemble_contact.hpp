#ifndef UTOPIA_ASSEMBLE_CONTACT_HPP
#define UTOPIA_ASSEMBLE_CONTACT_HPP 

#include <vector>
#include <utility>
#include <memory>


#include "utopia.hpp"
#include "libmesh/libmesh_common.h"

//forward decl
namespace moonolith {
	class Communicator;
}

namespace libMesh {
	class MeshBase;
	class DofMap;
}

namespace utopia {

	bool assemble_contact(
		moonolith::Communicator &comm,
		const std::shared_ptr<libMesh::MeshBase> &mesh,
		const std::shared_ptr<libMesh::DofMap> &dof_map,
		const unsigned int var_num,
		DSMatrixd &B,
		DSMatrixd &orthogonal_trafos,
		DVectord &gap,
		DSMatrixd &normals,
		DVectord &is_contact_node,
		const libMesh::Real search_radius,
		const std::vector< std::pair<int, int> > &tags,
		const bool use_biorth = true,
		const bool use_volume_differential = false);

	bool assemble_contact(
		moonolith::Communicator &comm,
		const std::shared_ptr<libMesh::MeshBase> &mesh,
		const std::shared_ptr<libMesh::DofMap> &dof_map,
		const unsigned int var_num,
		DSMatrixd &B,
		DSMatrixd &orthogonal_trafos,
		DVectord &gap,
		DSMatrixd &normals,
		DVectord &is_contact_node,
		const libMesh::Real search_radius,
		const int tag_1, 
		const int tag_2,
		const bool use_biorth = true,
		const bool use_volume_differential = false);


	bool assemble_contact(
		moonolith::Communicator &comm,
		const std::shared_ptr<libMesh::MeshBase> &mesh,
		const std::shared_ptr<libMesh::DofMap> &dof_map,
		const unsigned int var_num,
		DSMatrixd &B,
		DSMatrixd &orthogonal_trafos,
		DVectord &gap,
		DVectord &normals,
		DVectord &is_contact_node,
		const libMesh::Real search_radius,
		const std::vector< std::pair<int, int> > &tags,
		const bool use_biorth = true,
		const bool use_volume_differential = false);


	void convert_normal_matrix_to_vector(const DSMatrixd &mat, DVectord &vec);

}

#endif //UTOPIA_ASSEMBLE_CONTACT_HPP
