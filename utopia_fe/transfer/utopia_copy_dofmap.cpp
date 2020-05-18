
#include "utopia_copy_dofmap.hpp"

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

    Copy_Dof_Map::Copy_Dof_Map() {}

    void Copy_Dof_Map::copy_global_dofs(MeshBase &mesh,
                                        const std::shared_ptr<libMesh::DofMap> &original_dof_map,
                                        const unsigned int &var_num,
                                        std::vector<ElementDofMap> &dof_map,
                                        std::vector<ElementDofMap> &variable_type,
                                        std::vector<libMesh::dof_id_type> &handle_to_element_id) {
        const libMesh::dof_id_type n_elements = mesh.n_active_local_elem();

        std::vector<dof_id_type> temp;
        dof_map.resize(n_elements);
        handle_to_element_id.resize(n_elements);

        MeshBase::const_element_iterator e_it = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator e_end = mesh.active_local_elements_end();

        variable_type.resize(1);

        bool first = true;

        libMesh::dof_id_type local_element_id = 0;
        for (; e_it != e_end; ++e_it, ++local_element_id) {
            Elem *elem = *e_it;

            handle_to_element_id[local_element_id] = elem->id();

            original_dof_map->dof_indices(elem, temp, var_num);

            dof_map[local_element_id].global.insert(dof_map[local_element_id].global.end(), temp.begin(), temp.end());

            if (first) {
                variable_type[0].global.push_back(elem->type());
                first = false;
            }
        }
    }

    void Copy_Dof_Map::copy_var_order(DofMap &dofmap, std::vector<ElementDofMap> &variable_order) {
        variable_order.resize(1);
        FEType fe_type = dofmap.variable(0).type();
        variable_order[0].global.push_back(fe_type.order);
    }
}  // namespace utopia