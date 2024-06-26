#include "utopia_FESpacesAdapter.hpp"

#include "libmesh/elem.h"

#include <cmath>
#include <numeric>

using libMesh::dof_id_type;
using libMesh::DofMap;
using libMesh::Elem;
using libMesh::FEType;
using libMesh::MeshBase;
using libMesh::Node;

namespace utopia {

    FESpacesAdapter::FESpacesAdapter(const std::shared_ptr<MeshBase> &master,
                                     const std::shared_ptr<MeshBase> &slave,
                                     const std::shared_ptr<libMesh::DofMap> &dof_map_master,
                                     const std::shared_ptr<libMesh::DofMap> &dof_map_slave,
                                     const unsigned int &from_var_num,
                                     const unsigned int &to_var_num) {
        spaces_.reserve(2);
        spaces_.push_back(master);
        spaces_.push_back(slave);

        Copy_Dof_Map::copy_global_dofs(
            *master, dof_map_master, from_var_num, dof_maps_[0], var_type_[0], handle_to_element_id_[0]);
        Copy_Dof_Map::copy_global_dofs(
            *slave, dof_map_slave, to_var_num, dof_maps_[1], var_type_[1], handle_to_element_id_[1]);

        Copy_Dof_Map::copy_var_order(*dof_map_master, var_order_[0]);
        Copy_Dof_Map::copy_var_order(*dof_map_slave, var_order_[1]);
    }

}  // namespace utopia
