#include "utopia_LibMeshDofMapAdapterNew.hpp"


#define CHRONO_START() utopia::Chrono macro_chrono_; macro_chrono_.start();
#define CHRONO_END(macro_name_) { macro_chrono_.stop(); std::cout << macro_name_ << ":"<< macro_chrono_ << std::endl;}

namespace utopia {

    void ElementDofMapAdapter::init(
        libMesh::MeshBase &vol_mesh,
        const libMesh::DofMap &volume_dof_map,
        libMesh::MeshBase &surf_mesh,
        const libMesh::DofMap &surf_dof_map,
        const int var_num,
        const int surf_var_num)
    {
        CHRONO_START();
        //I do not think we need anything but 0 at the moment
        unsigned int comp = 0;
        unsigned int surf_comp = 0;
        unsigned int sys_num   = volume_dof_map.sys_number();
        unsigned int surf_sys_num = surf_dof_map.sys_number();

        std::vector<libMesh::dof_id_type> dof_surf;
        std::vector<libMesh::dof_id_type> dof_vol;

        const auto n_elem = surf_mesh.n_active_local_elem();
        dof_map_.resize(n_elem);

        SizeType local_el_idx = -1;

        for(const auto &b_elem : surf_mesh.active_local_element_ptr_range())
        {
            ++local_el_idx;
            
            auto &dof_object = dof_map_.dof_object(local_el_idx);
            dof_object.global_idx = b_elem->id();
            dof_object.block      = b_elem->subdomain_id();

            const libMesh::Elem * v_elem = b_elem->interior_parent();

            surf_dof_map.dof_indices(b_elem, dof_surf);
            auto n_dofs_x_el = dof_surf.size();
            
            dof_object.dofs.resize(n_dofs_x_el);

            // loop through all nodes in each boundary element.
            for (unsigned int node = 0; node < b_elem->n_nodes(); node++) {

                // Node in boundary element.
                const libMesh::Node * b_node = b_elem->node_ptr(node);

                for (unsigned int node_id = 0; node_id < v_elem->n_nodes(); node_id++)
                {
                    // Nodes in interior_parent element.
                    const libMesh::Node * v_node = v_elem->node_ptr(node_id);

                    const auto v_dof = v_node->dof_number(
                        sys_num,
                        var_num,
                        comp
                    );

                    if(v_node->absolute_fuzzy_equals(*b_node, 1e-14))
                    {
                        // Global dof_index for node in BoundaryMesh
                        const auto b_dof = b_node->dof_number(
                            surf_sys_num,
                            surf_var_num,
                            surf_comp);

                        auto dof_it = std::find(dof_surf.begin(), dof_surf.end(), b_dof);
                        assert(dof_it != dof_surf.end());

                        auto local_ind = std::distance(dof_surf.begin(), dof_it);
                        dof_object.dofs[local_ind] = v_dof;
                    }
                }
            }
        }

        CHRONO_END("boundary_element_node_permutation_map");
    }

}
