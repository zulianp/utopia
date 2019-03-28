#ifndef UTOPIA_COPY_DOFMAP_HPP
#define UTOPIA_COPY_DOFMAP_HPP

#include "utopia_ElementDofMap.hpp"

#include "libmesh/serial_mesh.h"
#include "libmesh/dof_map.h"

#include <vector>

namespace libMesh {
    class MeshBase;
    class DofMap;
}

namespace utopia {

 class Copy_Dof_Map{

    public:
     Copy_Dof_Map();

     static void copy_global_dofs(libMesh::MeshBase &space,
                                    const std::shared_ptr<libMesh::DofMap>  &original_dof_map,
                                    const unsigned int  &var_num,
                                    std::vector<ElementDofMap> &dof_map,
                                    std::vector<ElementDofMap> &variable_type,
                                    std::vector<libMesh::dof_id_type> &handle_two_element_id);

    static void copy_var_order(libMesh::DofMap &dofmap, std::vector<ElementDofMap> &variable_order);
    };
}

#endif //UTOPIA_COPY_DOFMAP_HPP