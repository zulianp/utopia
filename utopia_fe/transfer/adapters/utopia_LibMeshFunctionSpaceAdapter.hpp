#ifndef UTOPIA_LIBMESH_FUNCTION_SPACE_ADAPTER_HPP
#define UTOPIA_LIBMESH_FUNCTION_SPACE_ADAPTER_HPP

#include "utopia_LibMeshCollectionManager.hpp"
#include "utopia_LibMeshDofMapAdapter.hpp"

namespace utopia {

    class LibMeshFunctionSpaceAdapter {
    public:
        using ElementIter = typename libMesh::MeshBase::const_element_iterator;
        using Integer = libMesh::dof_id_type;

        class FEType {
        public:
            int family;
            int order;
        };

        LibMeshFunctionSpaceAdapter() : is_extracted_surface_(false) {}

        inline libMesh::MeshBase &mesh() {
            assert(mesh_);
            return *mesh_;
        }

        inline const libMesh::MeshBase &mesh() const {
            assert(mesh_);
            return *mesh_;
        }

        inline const moonolith::Storage<moonolith::Dofs> &element_dof_map() const { return dof_map_.element_dof_map(); }

        inline moonolith::Storage<moonolith::Dofs> &element_dof_map() { return dof_map_.element_dof_map(); }

        inline std::vector<libMesh::dof_id_type> &handle_to_element_id() { return handle_to_element_id_; }

        inline const std::vector<libMesh::dof_id_type> &handle_to_element_id() const { return handle_to_element_id_; }

        inline libMesh::dof_id_type handle_to_element_id(const std::size_t local_idx) const {
            assert(local_idx < handle_to_element_id_.size());
            return handle_to_element_id_[local_idx];
        }

        void set_mesh(const std::shared_ptr<libMesh::MeshBase> &mesh) { mesh_ = mesh; }

        void set_n_vars(const std::size_t n_vars) { fe_type_.resize(n_vars); }

        inline std::size_t n_vars() const { return fe_type_.size(); }

        inline int spatial_dim() const { return mesh_->spatial_dimension(); }

        FEType &fe_type(const std::size_t i) {
            assert(i < fe_type_.size());
            return fe_type_[i];
        }

        const FEType &fe_type(const std::size_t i) const {
            assert(i < fe_type_.size());
            return fe_type_[i];
        }

        const libMesh::FEType libmesh_fe_type(const std::size_t i) const {
            assert(i < fe_type_.size());
            const auto &f = fe_type(i);

            return libMesh::FEType(libMesh::Order(f.order), libMesh::FEFamily(f.family));
        }

        void boundary_ids_workaround(const libMesh::MeshBase &parent_mesh) {
            auto e_it = mesh_->active_local_elements_begin();
            auto e_end = mesh_->active_local_elements_end();

            for (; e_it != e_end; ++e_it) {
                auto *elem = *e_it;
                auto *parent = elem->interior_parent();
                assert(parent);

                std::size_t n_sides = parent->n_sides();

                auto c_elem = elem->centroid();

                bool found_side = false;
                for (std::size_t side_num = 0; side_num < n_sides; ++side_num) {
                    if (parent->neighbor_ptr(side_num) == nullptr) {
                        auto side_ptr = parent->build_side_ptr(side_num);
                        auto c_side = side_ptr->centroid();

                        if ((c_side - c_elem).norm() < 1e-14) {
                            // same element overwrite useless information
                            elem->subdomain_id() =
                                utopia::boundary_id(parent_mesh.get_boundary_info(), parent, side_num);
                            found_side = true;
                            break;
                        }
                    }
                }

                assert(found_side);
            }
        }

        void extract_surface_init_for_contact(const std::shared_ptr<libMesh::MeshBase> &mesh,
                                              const libMesh::DofMap &dof_map,
                                              const int var_num) {
            Chrono c;
            c.start();

            auto b_mesh = std::make_shared<libMesh::BoundaryMesh>(mesh->comm(), mesh->mesh_dimension() - 1);
            source_dof_map_ = utopia::make_ref(dof_map);

            mesh->get_boundary_info().sync(*b_mesh);

            es = utopia::make_unique<libMesh::EquationSystems>(*b_mesh);
            es->add_system<libMesh::LinearImplicitSystem>("boundary_sys");

            libMesh::FEType fe_type = dof_map.variable_type(var_num);
            auto &sys = es->get_system("boundary_sys");
            auto b_var_num = sys.add_variable("lambda", fe_type);
            surf_dof_map_ = make_ref(sys.get_dof_map());
            es->init();

            c.stop();

            std::cout << "libmesh get_boundary_info + eq system: " << c << std::endl;

            init_aux(b_mesh, sys.get_dof_map(), b_var_num);
            is_extracted_surface_ = true;
            boundary_ids_workaround(*mesh);

            // surf_mesh is extracted from vol_mesh
            dof_map_.init_for_contact(*mesh, dof_map, *b_mesh, sys.get_dof_map(), var_num, b_var_num);

            dof_map_.describe();
        }

        void init(const std::shared_ptr<libMesh::MeshBase> &mesh, const libMesh::DofMap &dof_map, const int var_num) {
            Chrono c;
            c.start();

            source_dof_map_ = utopia::make_ref(dof_map);
            libMesh::FEType fe_type = dof_map.variable_type(var_num);

            init_aux(mesh, dof_map, var_num);
            is_extracted_surface_ = false;

            dof_map_.init(*mesh, dof_map, var_num);

            dof_map_.describe();

            c.stop();
            std::cout << "LibMeshFunctionSpaceAdapter::init: " << c << std::endl;
        }

        // void extract_surface_init(
        //     const std::shared_ptr<libMesh::MeshBase> &mesh,
        //     const libMesh::DofMap &dof_map,
        //     const int var_num)
        // {

        //     Chrono c;
        //     c.start();

        //     auto b_mesh = std::make_shared<libMesh::BoundaryMesh>(mesh->comm(), mesh->mesh_dimension() - 1);
        //     source_dof_map_ = utopia::make_ref(dof_map);

        //     mesh->get_boundary_info().sync(*b_mesh);

        //     es = utopia::make_unique<libMesh::EquationSystems>(*b_mesh);
        //     es->add_system<libMesh::LinearImplicitSystem> ("boundary_sys");

        //     libMesh::FEType fe_type = dof_map.variable_type(var_num);
        //     auto &sys = es->get_system("boundary_sys");
        //     auto b_var_num = sys.add_variable("lambda", fe_type);
        //     surf_dof_map_ = make_ref(sys.get_dof_map());
        //     es->init();

        //     c.stop();

        //     std::cout << "libmesh get_boundary_info + eq system: " << c << std::endl;

        //     init_aux(b_mesh, sys.get_dof_map(), b_var_num);
        //     is_extracted_surface_ = true;
        //     boundary_ids_workaround(*mesh);

        //     //surf_mesh is extracted from vol_mesh
        //     dof_map_.init_for_surface(
        //         *mesh,
        //         dof_map,
        //         *b_mesh,
        //         sys.get_dof_map(),
        //         var_num,
        //         b_var_num);

        //     dof_map_.describe();
        // }

        static Integer tag(const ElementIter &e_it) { return (*e_it)->subdomain_id(); }

        void print_tags() const {
            auto e_it = mesh_->active_local_elements_begin();
            auto e_end = mesh_->active_local_elements_end();

            std::set<int> unique_tags;

            libMesh::dof_id_type local_element_id = 0;
            for (; e_it != e_end; ++e_it, ++local_element_id) {
                auto *elem = *e_it;
                unique_tags.insert(tag(e_it));
            }

            std::cout << "----------\n";
            std::cout << "tags:\n";

            for (auto t : unique_tags) {
                std::cout << t << " ";
            }

            std::cout << "----------\n";
        }

        const libMesh::Parallel::Communicator &comm() const { return this->mesh_->comm(); }

        const libMesh::DofMap &source_dof_map() const {
            assert(source_dof_map_);
            return *source_dof_map_;
        }

        void init_aux(const std::shared_ptr<libMesh::MeshBase> &mesh,
                      const libMesh::DofMap &dof_map,
                      const int var_num) {
            const libMesh::dof_id_type n_elements = mesh->n_active_local_elem();

            handle_to_element_id_.resize(n_elements);

            auto e_it = mesh->active_local_elements_begin();
            auto e_end = mesh->active_local_elements_end();

            libMesh::dof_id_type local_element_id = 0;
            for (; e_it != e_end; ++e_it, ++local_element_id) {
                auto *elem = *e_it;

                handle_to_element_id_[local_element_id] = elem->id();
            }

            // std::size_t n_vars = dof_map.n_variables();
            // this->set_n_vars(n_vars);

            // for(std::size_t i = 0; i < n_vars; ++i) {
            //     this->fe_type(i).family = dof_map.variable(i).type().family;
            //     this->fe_type(i).order  = dof_map.variable(i).type().order;
            // }

            // std::size_t n_vars = dof_map.n_variables();
            this->set_n_vars(1);

            this->fe_type(0).family = dof_map.variable(var_num).type().family;
            this->fe_type(0).order = dof_map.variable(var_num).type().order;

            mesh_ = mesh;
            is_extracted_surface_ = false;
        }

        // void init(
        //     const std::shared_ptr<libMesh::MeshBase> &mesh,
        //     const libMesh::DofMap &dof_map,
        //     const int var_num)
        // {
        //     const libMesh::dof_id_type n_elements = mesh->n_active_local_elem();

        //     std::vector<libMesh::dof_id_type> temp;
        //     dof_map_.resize(n_elements);
        //     handle_to_element_id_.resize(n_elements);

        //     auto e_it  = mesh->active_local_elements_begin();
        //     auto e_end = mesh->active_local_elements_end();

        //     libMesh::dof_id_type local_element_id = 0;
        //     for (; e_it != e_end; ++e_it, ++local_element_id){
        //         auto *elem = *e_it;

        //         handle_to_element_id_[local_element_id] = elem->id();

        //         dof_map.dof_indices(elem, temp, var_num);

        //         dof_map_[local_element_id].global_id = elem->id();
        //         dof_map_[local_element_id].global.insert(
        //             dof_map_[local_element_id].global.end(),
        //             temp.begin(),
        //             temp.end()
        //         );
        //     }

        //     // std::size_t n_vars = dof_map.n_variables();
        //     // this->set_n_vars(n_vars);

        //     // for(std::size_t i = 0; i < n_vars; ++i) {
        //     //     this->fe_type(i).family = dof_map.variable(i).type().family;
        //     //     this->fe_type(i).order  = dof_map.variable(i).type().order;
        //     // }

        //     std::size_t n_vars = dof_map.n_variables();
        //     this->set_n_vars(1);

        //     this->fe_type(0).family = dof_map.variable(var_num).type().family;
        //     this->fe_type(0).order  = dof_map.variable(var_num).type().order;

        //     mesh_ = mesh;
        //     is_extracted_surface_ = false;
        //     n_local_dofs_ = dof_map.n_local_dofs();
        // }

        const std::shared_ptr<const USparseMatrix> permutation() const { return dof_map_.permutation(); }

        const std::shared_ptr<const USparseMatrix> vector_permutation() const { return dof_map_.vector_permutation(); }

        // inline const moonolith::Storage<moonolith::Dofs> &element_dof_map() const
        // {
        //     return dof_map_.element_dof_map();
        // }

        // inline moonolith::Storage<moonolith::Dofs> &element_dof_map()
        // {
        //     return dof_map_.element_dof_map();
        // }

        // void bundary_permutation_matrix(
        //     const libMesh::MeshBase &boundary_mesh,
        //     const libMesh::DofMap &volume_dof_map,
        //     const libMesh::DofMap &surface_dof_map,
        //     unsigned int var_num,
        //     unsigned int comp,
        //     unsigned int b_var_num,
        //     unsigned int b_comp,
        //     USparseMatrix &matrix)
        // {
        //     bundary_permutation_map(
        //         boundary_mesh,
        //         volume_dof_map,
        //         surface_dof_map,
        //         var_num,
        //         comp,
        //         b_var_num,
        //         b_comp,
        //         boundary_to_volume_map
        //     );

        //     bundary_permutation_matrix_from_map(
        //         boundary_to_volume_map,
        //         *permutation_
        //     );
        // }

        // class Map {
        // public:
        //     SizeType from_range_begin, from_range_end;
        //     SizeType to_range_begin, to_range_end;
        //     std::vector<SizeType> idx;

        //     inline SizeType from_extent() const { return from_range_end - from_range_begin; }
        //     inline SizeType to_extent()   const { return to_range_end   - to_range_begin; }
        // };

        // static void bundary_permutation_matrix_from_map(
        //     const Map &map,
        //     USparseMatrix &mat)
        // {
        //     auto n_local_dof_vol = map.to_range_end - map.to_range_begin;
        //     auto n_local_dofs_surf = map.idx.size();

        //     mat = local_sparse(n_local_dof_vol, n_local_dofs_surf, 1);

        //     Write<USparseMatrix> w(mat);
        //     for(std::size_t i = 0; i < n_local_dofs_surf; ++i) {
        //         mat.set(map.idx[i], i + map.from_range_begin, 1.);
        //     }
        // }

        // static SizeType element_node_dof_map_and_permutation(
        //     const libMesh::MeshBase &mesh,
        //     const libMesh::DofMap &dof_map,
        //     const Map &map,
        //     std::vector<ElementDofMap> &elem_dof_map,
        //     USparseMatrix &mat)
        // {
        //     moonolith::Communicator comm(mesh.comm().get());

        //     auto n_local_dof_vol  = map.to_range_end   - map.to_range_begin;
        //     auto n_local_dof_surf = map.from_range_end - map.from_range_begin;

        //     SizeType dof_x_elem = 0, n_local_elems = mesh.n_active_local_elem();
        //     std::vector<libMesh::dof_id_type> dof_indices;

        //     {
        //         auto e_it = elements_begin(mesh);
        //         if(e_it != elements_end(mesh)) {
        //             dof_map.dof_indices(*e_it, dof_indices);
        //             dof_x_elem = dof_indices.size();
        //         }
        //     }

        //     std::vector<SizeType> dof_offsets(comm.size() + 1, 0);
        //     dof_offsets[comm.rank() + 1] = dof_x_elem * n_local_elems;

        //     comm.all_reduce(&dof_offsets[0], dof_offsets.size(), moonolith::MPISum());

        //     auto start_dof = dof_offsets[comm.rank()];
        //     auto n_local_e2n_dofs_surf = dof_offsets[comm.rank() + 1] - start_dof;

        //     mat = local_sparse(n_local_dof_vol, n_local_e2n_dofs_surf, 1);

        //     auto r_begin = map.from_range_begin;

        //     elem_dof_map.resize(n_local_elems);

        //     SizeType el_idx = 0;
        //     SizeType idx = start_dof;
        //     Write<USparseMatrix> w(mat);
        //     for(auto e_it = elements_begin(mesh); e_it != elements_end(mesh); ++e_it, ++el_idx) {
        //         dof_map.dof_indices(*e_it, dof_indices);

        //         auto &el_dof = elem_dof_map[el_idx];
        //         auto n = dof_indices.size();
        //         el_dof.global.resize(n);
        //         el_dof.global_id = (*e_it)->id();

        //         for(std::size_t k = 0; k < n; ++k, ++idx) {
        //             const auto i = dof_indices[k];
        //             const auto local_i = i - r_begin;
        //             assert(local_i < map.idx.size());

        //             mat.set(map.idx[local_i], idx, 1.);
        //             el_dof.global[k] = idx;
        //         }
        //     }

        //     return n_local_e2n_dofs_surf;
        // }

        // static void bundary_permutation_map(
        //     const libMesh::MeshBase &boundary_mesh,
        //     const libMesh::DofMap &volume_dof_map,
        //     const libMesh::DofMap &surface_dof_map,
        //     unsigned int var_num,
        //     unsigned int comp,
        //     unsigned int b_var_num,
        //     unsigned int b_comp,
        //     Map &map)
        // {
        //     std::vector<libMesh::dof_id_type> dof_boundary;

        //     unsigned int sys_num   = volume_dof_map.sys_number();
        //     unsigned int b_sys_num = surface_dof_map.sys_number();

        //     map.idx.resize(surface_dof_map.n_local_dofs());
        //     std::fill(map.idx.begin(), map.idx.end(), 0);

        //     Range rr(surface_dof_map.first_dof(), surface_dof_map.last_dof() + 1);

        //     map.from_range_begin = rr.begin();
        //     map.from_range_end = rr.end();

        //     map.to_range_begin = volume_dof_map.first_dof();
        //     map.to_range_end   = volume_dof_map.last_dof() + 1;

        //     // loop through all boundary elements.
        //     for (const auto & b_elem : boundary_mesh.active_local_element_ptr_range())
        //     {
        //         const libMesh::Elem * v_elem = b_elem->interior_parent();

        //         // loop through all nodes in each boundary element.
        //         for (unsigned int node = 0; node < b_elem->n_nodes(); node++) {

        //             // Node in boundary element.
        //             const libMesh::Node * b_node = b_elem->node_ptr(node);

        //             for (unsigned int node_id=0; node_id < v_elem->n_nodes(); node_id++)
        //             {
        //                 // Nodes in interior_parent element.
        //                 const libMesh::Node * v_node = v_elem->node_ptr(node_id);

        //                 const auto v_dof = v_node->dof_number(
        //                                                 sys_num,
        //                                                 var_num,
        //                                                 comp);

        //                 if (v_node->absolute_fuzzy_equals(*b_node, 1e-14))
        //                 {
        //                     // Global dof_index for node in BoundaryMesh
        //                     const auto b_dof = b_node->dof_number(
        //                                                     b_sys_num,
        //                                                     b_var_num,
        //                                                     b_comp);

        //                     if(rr.inside(b_dof)){
        //                         map.idx[b_dof - rr.begin()] = v_dof;
        //                     }
        //                 }
        //             }
        //         }
        //     }
        // }

        inline libMesh::dof_id_type n_local_elems() const { return mesh().n_active_local_elem(); }

        inline libMesh::dof_id_type n_local_dofs() const {
            // return n_local_dofs_;
            return dof_map_.n_local_dofs();
        }

        inline const libMesh::DofMap &surf_dof_map() const {
            assert(surf_dof_map_);
            return *surf_dof_map_;
        }

    private:
        std::shared_ptr<libMesh::MeshBase> mesh_;
        std::vector<libMesh::dof_id_type> handle_to_element_id_;
        std::vector<FEType> fe_type_;
        std::unique_ptr<libMesh::EquationSystems> es;

        LibMeshDofMapAdapter dof_map_;
        bool is_extracted_surface_;

        std::shared_ptr<const libMesh::DofMap> source_dof_map_;
        std::shared_ptr<const libMesh::DofMap> surf_dof_map_;
    };

    // class PermutationBuilder {
    // public:
    //     //1) build permutation from surface to volume
    //     //2) build permutation from element-node dofmap to volume
    //     //3) build permutation from surface element-node dofmap to volume

    template <class Bound, class Collection>
    class LibMeshElemenAdapter : public moonolith::ElementAdapterBase<Bound, Collection> {
    public:
        using super = moonolith::ElementAdapterBase<Bound, Collection>;
        using super::super;

        // inline const ElementDofMap &dofs() const
        // {
        //    assert(dofs_);
        //    return *dofs_;
        // }

        // void set_dofs(const ElementDofMap * dofs)
        // {
        //     dofs_ = dofs;
        // }

        inline const moonolith::Dofs &dofs() const {
            assert(dofs_);
            return *dofs_;
        }

        void set_dofs(const moonolith::Dofs *dofs) { dofs_ = dofs; }

    private:
        // const ElementDofMap * dofs_ = nullptr;
        const moonolith::Dofs *dofs_ = nullptr;
    };

}  // namespace utopia

namespace moonolith {

    template <class Bound>
    class ElementAdapter<Bound, utopia::LibMeshFunctionSpaceAdapter>
        : public utopia::LibMeshElemenAdapter<Bound, utopia::LibMeshFunctionSpaceAdapter> {
    public:
        using super = utopia::LibMeshElemenAdapter<Bound, utopia::LibMeshFunctionSpaceAdapter>;
        using super::super;
    };

    template <>
    class CollectionManager<utopia::LibMeshFunctionSpaceAdapter>
        : public utopia::LibMeshCollectionManager<utopia::LibMeshFunctionSpaceAdapter> {
    public:
        using super = utopia::LibMeshCollectionManager<utopia::LibMeshFunctionSpaceAdapter>;
        using super::super;
    };

}  // namespace moonolith

namespace utopia {
    using LibMeshCollectionManagerT = moonolith::CollectionManager<utopia::LibMeshFunctionSpaceAdapter>;
}

#endif  // UTOPIA_LIBMESH_FUNCTION_SPACE_ADAPTER_HPP
