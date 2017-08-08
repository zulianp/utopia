#ifndef UTOPIA_NORMAL_TANGENTIAL_COORDINATE_SYSTEM_HPP
#define UTOPIA_NORMAL_TANGENTIAL_COORDINATE_SYSTEM_HPP

#include "utopia.hpp"
#include <vector>

namespace libMesh {
	class MeshBase;
	class DofMap;
}

namespace utopia {

	void plot_scaled_normal_field(
		const libMesh::MeshBase &mesh,
		const DVectord &normals,
		const DVectord &scale,
		const std::string &name = "normal_field");

	bool assemble_normal_tangential_transformation(
		const libMesh::MeshBase &mesh, 
		const libMesh::DofMap &dof_map,
		const std::vector<int> &boundary_tags, 
		DVectord &is_normal_component,
		DVectord &normals,
		DSMatrixd &mat);
}

#endif //UTOPIA_NORMAL_TANGENTIAL_COORDINATE_SYSTEM_HPP
