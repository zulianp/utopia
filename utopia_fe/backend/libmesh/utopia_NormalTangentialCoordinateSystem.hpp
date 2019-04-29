#ifndef UTOPIA_NORMAL_TANGENTIAL_COORDINATE_SYSTEM_HPP
#define UTOPIA_NORMAL_TANGENTIAL_COORDINATE_SYSTEM_HPP

#include "utopia_fe_base.hpp"
#include <vector>

namespace libMesh {
    class MeshBase;
    class DofMap;
}

namespace utopia {

    void scale_normal_vector_with_gap(const int dim, const UVector &normals, const UVector &gap, UVector &out);

    bool assemble_normal_tangential_transformation(
        const libMesh::MeshBase &mesh,
        const libMesh::DofMap &dof_map,
        const std::vector<int> &boundary_tags,
        UVector &is_normal_component,
        UVector &normals,
        USparseMatrix &mat);
}

#endif //UTOPIA_NORMAL_TANGENTIAL_COORDINATE_SYSTEM_HPP
