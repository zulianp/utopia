#include "utopia_LibMeshDofMapAdapter.hpp"

namespace utopia {

#define CHRONO_START() utopia::Chrono macro_chrono_; macro_chrono_.start();
#define CHRONO_END(macro_name_) { macro_chrono_.stop(); std::cout << macro_name_ << ":"<< macro_chrono_ << std::endl;}


    static std::size_t max_nnz_x_row(const libMesh::DofMap &dof_map)
    {
        std::size_t nnz = 0;

        if(!dof_map.get_n_nz().empty()) {
            nnz = *std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end());

        }

        if(!dof_map.get_n_oz().empty()) {
            nnz += *std::max_element(dof_map.get_n_oz().begin(), dof_map.get_n_oz().end());
        }

        return nnz;
    }

    //surf_mesh is extracted from vol_mesh
    void LibMeshDofMapAdapter::init_for_contact(
        libMesh::MeshBase &vol_mesh,
        const libMesh::DofMap &vol_dof_map,
        libMesh::MeshBase &surf_mesh,
        const libMesh::DofMap &surf_dof_map,
        const int var_num,
        const int surf_var_num)
    {
        //scalar dofs
        permutation_        = std::make_shared<USparseMatrix>();
        
        //vector dofs
        vector_permutation_ = std::make_shared<USparseMatrix>();

        SizeType spatial_dim = surf_mesh.spatial_dimension();
        max_nnz_ = max_nnz_x_row(surf_dof_map);

        boundary_element_node_permutation_map(
            surf_mesh,
            vol_dof_map,
            surf_dof_map,
            var_num,
            surf_var_num,
            mapping_
        );

        element_node_dof_map_and_permutation(
                               surf_mesh,
                               surf_dof_map,
                               mapping_,
                               dofs_,
                                *permutation_
                            );

        //vector permuation
        vector_permuation_map_from_map(spatial_dim, mapping_, vector_mapping_);

        permutation_matrix_from_map(
                        vector_mapping_,
                        *vector_permutation_
                    );
    }

    void LibMeshDofMapAdapter::init_for_surface(
        libMesh::MeshBase &vol_mesh,
        const libMesh::DofMap &vol_dof_map,
        libMesh::MeshBase &surf_mesh,
        const libMesh::DofMap &surf_dof_map,
        const int var_num,
        const int surf_var_num)
    {
        //scalar dofs
        permutation_        = std::make_shared<USparseMatrix>();
        
        //vector dofs
        // vector_permutation_ = std::make_shared<USparseMatrix>();

        SizeType spatial_dim = surf_mesh.spatial_dimension();
        max_nnz_ = max_nnz_x_row(surf_dof_map);

        boundary_permutation_map(
            surf_mesh,
            vol_dof_map,
            surf_dof_map,
            var_num,
            surf_var_num,
            mapping_
        );

        permutation_matrix_from_map(
                               mapping_,
                                *permutation_
                            );

        //vector permuation
        // vector_permuation_map_from_map(spatial_dim, mapping_, vector_mapping_);

        // permutation_matrix_from_map(
        //                 vector_mapping_,
        //                 *vector_permutation_
        //             );
    }

    void init(
        const std::shared_ptr<libMesh::MeshBase> &mesh,
        const libMesh::DofMap &dof_map,
        const int var_num,
        const bool build_element_dof_map)
    {
        assert(false);
    }

    //let us assume that the displacement degees of freedom are consecutive 
    //in the volume dofmap
    void LibMeshDofMapAdapter::vector_permuation_map_from_map(
        const std::size_t spatial_dim,
        const Mapping &mapping,
        Mapping &vector_mapping)
    {
        assert(spatial_dim > 0);

        const auto n_mappings = mapping.size();
        vector_mapping.idx.resize(n_mappings * spatial_dim);

        for(std::size_t i = 0; i < n_mappings; ++i) {
            const std::size_t i_offset = i * spatial_dim;
            
            for(std::size_t d = 0; d < spatial_dim; ++d) {
                assert(i_offset + d < n_mappings * spatial_dim);
                vector_mapping.idx[i_offset + d] = mapping.idx[i] + d;
                // std::cout << "map[" << (i_offset + d) <<  "] = " << vector_mapping.idx[i_offset + d] << std::endl;;
            }
        }

        vector_mapping.from_range_begin = mapping.from_range_begin * spatial_dim;
        vector_mapping.from_range_end   = mapping.from_range_end   * spatial_dim;

        vector_mapping.to_range_begin   = mapping.to_range_begin;
        vector_mapping.to_range_end     = mapping.to_range_end;
    }

    void LibMeshDofMapAdapter::permutation_matrix_from_map(
        const Mapping &map,
        USparseMatrix &mat)
    {
        CHRONO_START();

        auto n_local_dof_vol = map.to_extent();
        auto n_mappings      = map.idx.size();

        mat = local_sparse(n_local_dof_vol, map.from_extent(), max_nnz_);

        {
            Write<USparseMatrix> w(mat);

            for(std::size_t i = 0; i < n_mappings; ++i) {
                mat.set(map.idx[i], i + map.from_range_begin, 1.);
            }
        }

        CHRONO_END("permutation_matrix_from_map");
    }

    void LibMeshDofMapAdapter::boundary_permutation_map(
        const libMesh::MeshBase &surf_mesh,
        const libMesh::DofMap   &volume_dof_map,
        const libMesh::DofMap   &surf_dof_map,
        unsigned int var_num,
        unsigned int surf_var_num,
        Mapping &map)
    {
        CHRONO_START();

        //I do not think we need anything but 0 at the moment
        unsigned int comp = 0;
        unsigned int surf_comp = 0;

        std::vector<libMesh::dof_id_type> dof_boundary;

        unsigned int sys_num   = volume_dof_map.sys_number();
        unsigned int surf_sys_num = surf_dof_map.sys_number();

        map.idx.resize(surf_dof_map.n_local_dofs());
        std::fill(map.idx.begin(), map.idx.end(), 0);

        Range rr(surf_dof_map.first_dof(), surf_dof_map.last_dof() + 1);

        map.from_range_begin = rr.begin();
        map.from_range_end   = rr.end();

        map.to_range_begin = volume_dof_map.first_dof();
        map.to_range_end   = volume_dof_map.last_dof() + 1;

        // loop through all boundary elements.
        for (const auto &b_elem : surf_mesh.active_local_element_ptr_range())
        {
            const libMesh::Elem * v_elem = b_elem->interior_parent();

            // loop through all nodes in each boundary element.
            for (unsigned int node = 0; node < b_elem->n_nodes(); node++) {
                
                // Node in boundary element.
                const libMesh::Node * b_node = b_elem->node_ptr(node);

                for (unsigned int node_id=0; node_id < v_elem->n_nodes(); node_id++)
                {
                    // Nodes in interior_parent element.
                    const libMesh::Node * v_node = v_elem->node_ptr(node_id);

                    const auto v_dof = v_node->dof_number(
                                                    sys_num,
                                                    var_num,
                                                    comp);

                    if (v_node->absolute_fuzzy_equals(*b_node, 1e-14))
                    {
                        // Global dof_index for node in BoundaryMesh
                        const auto b_dof = b_node->dof_number(
                                                        surf_sys_num,
                                                        surf_var_num,
                                                        surf_comp);
                       
                        if(rr.inside(b_dof)){
                            map.idx[b_dof - rr.begin()] = v_dof;
                        }
                    }
                }
            }
        }

        CHRONO_END("boundary_permutation_map");
    }

    void LibMeshDofMapAdapter::boundary_element_node_permutation_map(
        const libMesh::MeshBase &surf_mesh,
        const libMesh::DofMap   &volume_dof_map,
        const libMesh::DofMap   &surf_dof_map,
        unsigned int var_num,
        unsigned int surf_var_num,
        Mapping &map)
    {
        CHRONO_START();
        //I do not think we need anything but 0 at the moment
        unsigned int comp = 0;
        unsigned int surf_comp = 0;

        std::vector<libMesh::dof_id_type> dof_boundary;

        unsigned int sys_num   = volume_dof_map.sys_number();
        unsigned int surf_sys_num = surf_dof_map.sys_number();

        Range rr = element_node_range(
                surf_mesh,
                surf_dof_map
        );

        map.idx.resize(rr.extent());
        std::fill(map.idx.begin(), map.idx.end(), 0);

        map.from_range_begin = rr.begin();
        map.from_range_end   = rr.end();

        map.to_range_begin = volume_dof_map.first_dof();
        map.to_range_end   = volume_dof_map.last_dof() + 1;

        std::vector<libMesh::dof_id_type> dof_indices;
        // loop through all boundary elements.

        SizeType local_el_idx = -1;
        for(const auto &b_elem : surf_mesh.active_local_element_ptr_range())
        {
            ++local_el_idx;
            const libMesh::Elem * v_elem = b_elem->interior_parent();

            surf_dof_map.dof_indices(b_elem, dof_indices);
            auto n_dofs_x_el = dof_indices.size();

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
                                                    comp);

                    if (v_node->absolute_fuzzy_equals(*b_node, 1e-14))
                    {
                        // Global dof_index for node in BoundaryMesh
                        const auto b_dof = b_node->dof_number(
                                                        surf_sys_num,
                                                        surf_var_num,
                                                        surf_comp);
                       
                        auto dof_it = std::find(dof_indices.begin(), dof_indices.end(), b_dof);
                        assert(dof_it != dof_indices.end());

                        auto local_ind = std::distance(dof_indices.begin(), dof_it);
                        auto map_ind = local_el_idx * n_dofs_x_el + local_ind;

                        assert(map.idx[map_ind] == 0);
                        map.idx[map_ind] = v_dof;
                    }
                }
            }
        }

        CHRONO_END("boundary_element_node_permutation_map");
    }

    SizeType LibMeshDofMapAdapter::count_dof_x_elem(
        const libMesh::MeshBase &mesh,
        const libMesh::DofMap &dof_map
        )
    {
        SizeType dof_x_elem = 0;
        std::vector<libMesh::dof_id_type> dof_indices;

        auto e_it = elements_begin(mesh);
        if(e_it != elements_end(mesh)) {
            dof_map.dof_indices(*e_it, dof_indices);
            dof_x_elem = dof_indices.size();
        }

        return dof_x_elem;
    }

    Range LibMeshDofMapAdapter::element_node_range(
        const libMesh::MeshBase &mesh,
        const libMesh::DofMap &dof_map)
    {
        moonolith::Communicator comm(mesh.comm().get());

        SizeType dof_x_elem = count_dof_x_elem(mesh, dof_map), n_local_elems = mesh.n_active_local_elem();
        std::vector<libMesh::dof_id_type> dof_indices;

        std::vector<SizeType> dof_offsets(comm.size() + 1, 0);
        dof_offsets[comm.rank() + 1] = dof_x_elem * n_local_elems;

        comm.all_reduce(&dof_offsets[0], dof_offsets.size(), moonolith::MPISum());

        return Range(
            dof_offsets[comm.rank()],
            dof_offsets[comm.rank() + 1]);
    }

    void LibMeshDofMapAdapter::element_node_dof_map_and_permutation(
        const libMesh::MeshBase &mesh,
        const libMesh::DofMap &dof_map,
        const Mapping &map,
        std::vector<ElementDofMap> &elem_dof_map,
        USparseMatrix &mat)
    {
        CHRONO_START();

        moonolith::Communicator comm(mesh.comm().get());

        auto n_local_dof_vol  = map.to_range_end - map.to_range_begin;            
        auto n_local_elems = mesh.n_active_local_elem();

        mat = local_sparse(n_local_dof_vol, map.from_extent(), max_nnz_);

        auto r_begin = map.from_range_begin;

        elem_dof_map.resize(n_local_elems);
        std::vector<libMesh::dof_id_type> dof_indices;

        SizeType el_idx = 0;
        SizeType idx = map.from_range_begin;
        Write<USparseMatrix> w(mat);
       
        for(auto e_it = elements_begin(mesh); e_it != elements_end(mesh); ++e_it, ++el_idx) {
            dof_map.dof_indices(*e_it, dof_indices);

            auto &el_dof = elem_dof_map[el_idx];
            auto n = dof_indices.size();
            el_dof.global.resize(n);
            el_dof.global_id = (*e_it)->id();

            for(std::size_t k = 0; k < n; ++k, ++idx) {
                const auto i = dof_indices[k];
                const auto local_i = i - r_begin;
                assert(local_i < map.idx.size());

                mat.set(map.idx[idx], idx, 1.);
                el_dof.global[k] = idx;
            }
        }

        CHRONO_END("element_node_dof_map_and_permutation");
    }


}
