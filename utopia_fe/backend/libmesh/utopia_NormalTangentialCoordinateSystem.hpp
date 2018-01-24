#ifndef UTOPIA_NORMAL_TANGENTIAL_COORDINATE_SYSTEM_HPP
#define UTOPIA_NORMAL_TANGENTIAL_COORDINATE_SYSTEM_HPP

#include "utopia.hpp"
#include <vector>

namespace libMesh {
	class MeshBase;
	class DofMap;
}

namespace utopia {
	
	void scale_normal_vector_with_gap(const int dim, const DVectord &normals, const DVectord &gap, DVectord &out);
	
	bool assemble_normal_tangential_transformation(
		const libMesh::MeshBase &mesh, 
		const libMesh::DofMap &dof_map,
		const std::vector<int> &boundary_tags, 
		DVectord &is_normal_component,
		DVectord &normals,
		DSMatrixd &mat);
}

#endif //UTOPIA_NORMAL_TANGENTIAL_COORDINATE_SYSTEM_HPP
