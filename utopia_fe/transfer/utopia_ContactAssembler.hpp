#ifndef UTOPIA_CONTACT_ASSEMBLER_HPP
#define UTOPIA_CONTACT_ASSEMBLER_HPP

#include "utopia_TensorProductWithIdentity.hpp"

#include "libmesh/mesh.h"
#include "libmesh/elem.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/boundary_mesh.h"
#include "libmesh/dof_map.h"

#include "moonolith_mesh_adapter.hpp"
#include "moonolith_polygon.hpp"
#include "moonolith_line.hpp"
#include "moonolith_check_stream.hpp"
#include "moonolith_make_unique.hpp"
#include "moonolith_polyhedron.hpp"
#include "moonolith_polygon.hpp"
#include "moonolith_line.hpp"
#include "moonolith_map_quadrature.hpp"
#include "moonolith_affine_transform.hpp"
#include "utopia_LibMeshToMoonolithConvertions.hpp"

#include "utopia_ElementDofMap.hpp"
#include "MortarAssemble.hpp"

#include <cassert>

namespace utopia {

    template<class FunctionSpace>
    class LibMeshCollectionManager {
    public:
        using Elem  = libMesh::Elem;
        using Point = libMesh::Point;
        using ElementIter = typename libMesh::MeshBase::const_element_iterator;
        using Integer = libMesh::dof_id_type;

        const libMesh::Parallel::Communicator &comm;

        LibMeshCollectionManager(const libMesh::Parallel::Communicator &comm)
        : comm(comm)
        {}

        //we need to add this information because of the lack of local indexing in libmesh meshes
        template<class Adapter>
        static void data_added(const Integer local_idx, Adapter &a)
        {
            const auto &space = a.collection();
            assert( space.handle_to_element_id(local_idx) == a.elem().id() );
            a.set_dofs(&space.dof_map()[local_idx]);
        }

        static const Elem &elem(const FunctionSpace &space, const ElementIter &e_it) 
        {
            return **e_it;
        }

        static const Elem &elem(const FunctionSpace &space, const Integer idx) 
        {
            return *space.mesh().elem(idx);
        }

        static Integer tag(const FunctionSpace &space, const ElementIter &e_it)
        {
            return space.tag(e_it);
        }

        static Integer n_elements(const FunctionSpace &space)
        {
            if(space.mesh().is_serial()) {
                auto n_elems = space.mesh().n_active_local_elem();
                if(n_elems == 0) {
                    return space.mesh().n_elem();
                } else {
                    return n_elems;
                }
            } else {
                return space.mesh().n_active_local_elem();
            }
        }

        static ElementIter elements_begin(const FunctionSpace &space)
        {
            if(space.mesh().is_serial()) {
                auto n_elems = space.mesh().n_active_local_elem();

                if(n_elems == 0) {
                    return space.mesh().elements_begin();
                } else {
                    return space.mesh().active_local_elements_begin();
                }
            }

            return space.mesh().active_local_elements_begin();
        }

        static ElementIter elements_end(const FunctionSpace &space)
        {
          if(space.mesh().is_serial()) {
              auto n_elems = space.mesh().n_active_local_elem();

              if(n_elems == 0) {
                  return space.mesh().elements_end();
              } else {
                  return space.mesh().active_local_elements_end();
              }
          }

          return space.mesh().active_local_elements_end();
        }

        static Integer handle(const FunctionSpace &space, const ElementIter &e_it)
        {
            return (*e_it)->id();
        }

        static bool skip(const FunctionSpace &space, const ElementIter &e_it)
        {
            return false;
        }

        template<class Bound>
        static void fill_bound(const FunctionSpace &space, const Integer handle, Bound &bound, const double blow_up)
        {
            const auto &mesh = space.mesh();
            const auto &e = *mesh.elem(handle);
            const auto dim = mesh.spatial_dimension();

            if(dim > mesh.mesh_dimension()) {
                std::vector<double> p(dim);
                
                Point nn;
                compute_side_normal(dim, e, nn);

                for(Integer i = 0; i < e.n_nodes(); ++i) {
                    
                    const auto &q = e.node_ref(i);

                    for(Integer d = 0; d < dim; ++d) {
                        p[d] = q(d) + blow_up * nn(d);
                    }

                    bound += p;

                    for(Integer d = 0; d < dim; ++d) {
                        p[d] = q(d) - blow_up * nn(d);
                    }

                    bound += p;
                }

            } else {

                std::vector<double> p(dim);
                
                for(Integer i = 0; i < e.n_nodes(); ++i) {
                    const auto &q = e.node_ref(i);
                    
                    for(Integer d = 0; d < dim; ++d) {
                        p[d] = q(d);
                    }

                    bound += p;
                }
            }
        }

        template<class Iterator>
        static void serialize(const FunctionSpace &space, const Iterator &begin, const Iterator &end, moonolith::OutputStream &os)
        {
            CHECK_STREAM_WRITE_BEGIN("serialize", os);

            const auto &mesh = space.mesh();

            const int dim         = mesh.spatial_dimension();
            const long n_elements = std::distance(begin, end);
            const auto &dof_map   = space.dof_map();
            // const auto &fe_types  = space.fe_types();

            std::set<long> nodeIds;
            std::map<long, long> mapping;

            for(Iterator it = begin; it != end; ++it) {

                const libMesh::dof_id_type local_element_id = *it;
                const libMesh::dof_id_type global_element_id = space.handle_to_element_id(local_element_id);

                const libMesh::Elem *elem = mesh.elem(global_element_id);

                for(libMesh::dof_id_type j = 0; j != elem->n_nodes(); ++j) {
                    nodeIds.insert(elem->node(j));
                }

            }

            long n_nodes = nodeIds.size();

            // Estimate for allocation
            os.request_space( (n_elements * 8 + n_nodes * dim) * (sizeof(double) + sizeof(long)) );

            //WRITE 1
            os << dim;// << role;


            int index = 0;
            for (auto nodeId : nodeIds) {
                mapping[nodeId] = index++;
            }

            //WRITE 2
            os << n_nodes;

            //WRITE 6
            os << n_elements;

            for(auto node_id : nodeIds){

                const libMesh::Point &p = mesh.node(node_id);

                for(int i = 0; i < LIBMESH_DIM; ++i) {

                    //WRITE 3
                    os << p(i);

                }

            }

            std::vector<libMesh::dof_id_type> indices_vector;

            CHECK_STREAM_WRITE_BEGIN("elements", os);

            for(Iterator it = begin; it != end; ++it) {
                const libMesh::dof_id_type local_element_id = *it;
                const libMesh::dof_id_type global_element_id = space.handle_to_element_id(local_element_id);

                const libMesh::Elem *elem = mesh.elem(global_element_id);

                const int e_n_nodes = elem->n_nodes();

                const int type = elem->type();

                //WRITE 7
                os << type << e_n_nodes;

                for (int i = 0; i != e_n_nodes; ++i) {

                    auto it = mapping.find(elem->node(i));

                    assert(it != mapping.end());

                    int index = it->second;

                    //WRITE 8
                    os << index;

                }

                //WRITE 9
                assert(!dof_map.at(local_element_id).empty());
                os << dof_map.at(local_element_id);
                //WRITE 10
                int volume_tag = elem->subdomain_id();
                os << volume_tag;
            }

            CHECK_STREAM_WRITE_END("elements", os);



            //WRITE 11
            int n_vars = space.n_vars();

            os << n_vars;
            for(int i = 0; i < n_vars; ++i) {
                int fe_family = space.fe_type(i).family;
                int fe_order  = space.fe_type(i).order;
                os << fe_family << fe_order;
            }


            CHECK_STREAM_WRITE_END("serialize", os);
        }

        std::unique_ptr<FunctionSpace> build(moonolith::InputStream &is) const
        {

            CHECK_STREAM_READ_BEGIN("serialize", is);

            auto space = moonolith::make_unique<FunctionSpace>();
            auto &dof_map = space->dof_map();
            auto &handle_to_element_id = space->handle_to_element_id();

            using namespace std;

            //READ 1
            int dim, role;
            is >> dim;// >> role;

            //READ 2
            long n_nodes;
            is >> n_nodes;

            //READ 6
            long n_elements;
            is >> n_elements;

            auto mesh_ptr = std::make_shared<libMesh::SerialMesh>(comm, dim);

            mesh_ptr->reserve_nodes(n_nodes);

            for (long iii = 0; iii < n_nodes; ++iii) {

                libMesh::Point p;

                for(int j = 0; j < LIBMESH_DIM; ++j) {
                    //READ 3
                    is >> p(j);
                }

                mesh_ptr->add_point(p);
            }

            dof_map.resize(n_elements);

            handle_to_element_id.resize(n_elements);


            CHECK_STREAM_READ_BEGIN("elements", is);

            for(long i = 0; i < n_elements; ++i) {
                handle_to_element_id[i] = i;
                //READ 7

                int type, e_n_nodes;

                is >> type >> e_n_nodes;

                auto elem = libMesh::Elem::build(libMesh::ElemType(type)).release();


                int index;

                for (int ii = 0; ii != e_n_nodes; ++ii) {

                    //READ 8
                    is >> index;

                    elem->set_node(ii) = & mesh_ptr->node(index);
                }

                mesh_ptr->add_elem(elem);

                libmesh_assert(elem);

                //READ 9
                is >> dof_map.at(i);

                //READ 10
                int volume_tag;

                is >> volume_tag;

                elem->subdomain_id()=volume_tag;
            }

            CHECK_STREAM_READ_END("elements", is);

            //READ 11
            int n_vars;

            is >> n_vars;
            space->set_n_vars(n_vars);

            for(int i = 0; i < n_vars; ++i) {
                int fe_family;
                int fe_order;
               
                is >> fe_family >> fe_order;

                space->fe_type(i).family = fe_family;
                space->fe_type(i).order = fe_order;
            }

            space->set_mesh(mesh_ptr);

            CHECK_STREAM_READ_END("serialize", is);
            return space;
        }

    };

    class LibMeshFunctionSpaceAdapter {
    public:

        using ElementIter = typename libMesh::MeshBase::const_element_iterator;
        using Integer = libMesh::dof_id_type;

        class FEType {
        public:
            int family;
            int order;
        };

        LibMeshFunctionSpaceAdapter() : is_extracted_surface_(false)
        {}

        inline libMesh::MeshBase &mesh()
        {
            assert(mesh_);
            return *mesh_;
        }

        inline const libMesh::MeshBase &mesh() const
        {
            assert(mesh_);
            return *mesh_;
        }

        inline std::vector<ElementDofMap> &dof_map() 
        {
            return dof_map_;
        }

        inline const std::vector<ElementDofMap> &dof_map() const
        {
            return dof_map_;
        }

        inline std::vector<libMesh::dof_id_type> &handle_to_element_id()
        {
            return handle_to_element_id_;
        }

        inline const std::vector<libMesh::dof_id_type> &handle_to_element_id() const
        {
            return handle_to_element_id_;
        }

        inline libMesh::dof_id_type handle_to_element_id(const std::size_t local_idx) const
        {
            assert(local_idx < handle_to_element_id_.size());
            return handle_to_element_id_[local_idx];
        }

        void set_mesh(const std::shared_ptr<libMesh::MeshBase> &mesh)
        {
            mesh_ = mesh;
        }

        void set_n_vars(const std::size_t n_vars)
        {
            fe_type_.resize(n_vars);
        }

        inline std::size_t n_vars() const 
        {
            return fe_type_.size();
        } 

        FEType &fe_type(const std::size_t i)
        {
            assert(i < fe_type_.size());
            return fe_type_[i];
        }

        const FEType &fe_type(const std::size_t i) const
        {
            assert(i < fe_type_.size());
            return fe_type_[i];
        }

        const libMesh::FEType libmesh_fe_type(const std::size_t i) const
        {
            assert(i < fe_type_.size());
            const auto &f = fe_type(i);

            return libMesh::FEType(
                libMesh::Order(f.order),
                libMesh::FEFamily(f.family)
            );
        }

        void boundary_ids_workaround(const libMesh::MeshBase &parent_mesh)
        {
            auto e_it  = mesh_->active_local_elements_begin();
            auto e_end = mesh_->active_local_elements_end();

            for (; e_it != e_end; ++e_it){
                auto *elem = *e_it;
                auto *parent = elem->interior_parent();
                assert(parent);

                std::size_t n_sides = parent->n_sides();

                auto c_elem = elem->centroid();

                bool found_side = false;
                for(std::size_t side_num = 0; side_num < n_sides; ++side_num) {
                    if(parent->neighbor_ptr(side_num) == nullptr) {
                        auto side_ptr = parent->build_side_ptr(side_num);
                        auto c_side   = side_ptr->centroid();

                        if((c_side - c_elem).norm() < 1e-14) {
                            //same element overwrite useless information
                            elem->subdomain_id() = parent_mesh.get_boundary_info().boundary_id(parent, side_num);
                            found_side = true;
                            break;
                        }
                    }
                }

                assert(found_side);
            }
        }

        void extract_surface_init(
            const std::shared_ptr<libMesh::MeshBase> &mesh,
            const libMesh::DofMap &dof_map,
            const int var_num)
        {
            auto b_mesh = std::make_shared<libMesh::BoundaryMesh>(mesh->comm(), mesh->mesh_dimension() - 1);

            mesh->get_boundary_info().sync(*b_mesh);

            es = utopia::make_unique<libMesh::EquationSystems>(*b_mesh);
            es->add_system<libMesh::LinearImplicitSystem> ("boundary_sys");

            libMesh::FEType fe_type = dof_map.variable_type(var_num);
            auto &sys = es->get_system("boundary_sys");
            auto b_var_num = sys.add_variable("lambda", fe_type); 
            es->init();

            init(b_mesh, sys.get_dof_map(), b_var_num);
            is_extracted_surface_ = true;
            boundary_ids_workaround(*mesh);

            permutation_ = std::make_shared<USparseMatrix>();

            // bundary_permutation_matrix(
            // *b_mesh,
            // dof_map,
            // sys.get_dof_map(),
            // var_num,
            // 0,
            // 0,
            // 0,
            // *permutation_);

            bundary_permutation_map(
                *b_mesh,
                dof_map,
                sys.get_dof_map(),
                var_num,
                0,
                0,
                0,
                boundary_to_volume_map
            );

            bundary_permutation_matrix_from_map(
                boundary_to_volume_map,
                *permutation_
            );

            n_local_dofs_ =  sys.get_dof_map().n_local_dofs();
        }

        static Integer tag(const ElementIter &e_it)
        {
            return (*e_it)->subdomain_id();
        }

        void print_tags() const
        {
            auto e_it  = mesh_->active_local_elements_begin();
            auto e_end = mesh_->active_local_elements_end();

            std::set<int> unique_tags;

            libMesh::dof_id_type local_element_id = 0;
            for (; e_it != e_end; ++e_it, ++local_element_id){
                auto *elem = *e_it;
                unique_tags.insert(tag(e_it));
            }

            std::cout << "----------\n";
            std::cout << "tags:\n";

            for(auto t : unique_tags) {
                std::cout << t << " ";
            }

            std::cout << "----------\n";

        }

        const libMesh::Parallel::Communicator &comm() const
        {
            return this->mesh_->comm(); 
        }

        void init(
            const std::shared_ptr<libMesh::MeshBase> &mesh,
            const libMesh::DofMap &dof_map,
            const int var_num)
        {
            const libMesh::dof_id_type n_elements = mesh->n_active_local_elem();

            std::vector<libMesh::dof_id_type> temp;
            dof_map_.resize(n_elements);
            handle_to_element_id_.resize(n_elements);

            auto e_it  = mesh->active_local_elements_begin();
            auto e_end = mesh->active_local_elements_end();

            libMesh::dof_id_type local_element_id = 0;
            for (; e_it != e_end; ++e_it, ++local_element_id){
                auto *elem = *e_it;

                handle_to_element_id_[local_element_id] = elem->id();

                dof_map.dof_indices(elem, temp, var_num);

                dof_map_[local_element_id].global_id = elem->id();
                dof_map_[local_element_id].global.insert(
                    dof_map_[local_element_id].global.end(),
                    temp.begin(),
                    temp.end()
                );
            }

            // std::size_t n_vars = dof_map.n_variables();
            // this->set_n_vars(n_vars);

            // for(std::size_t i = 0; i < n_vars; ++i) {
            //     this->fe_type(i).family = dof_map.variable(i).type().family;
            //     this->fe_type(i).order  = dof_map.variable(i).type().order;
            // }

            std::size_t n_vars = dof_map.n_variables();
            this->set_n_vars(1);

            this->fe_type(0).family = dof_map.variable(var_num).type().family;
            this->fe_type(0).order  = dof_map.variable(var_num).type().order;

            mesh_ = mesh;
            is_extracted_surface_ = false;
            n_local_dofs_ = dof_map.n_local_dofs();
        }

        const std::shared_ptr<USparseMatrix> &permutation()
        {
            return permutation_;
        }

        const std::shared_ptr<USparseMatrix> &tensor_permutation()
        {
            return tensor_permutation_;
        }

        static void bundary_permutation_matrix(
            const libMesh::MeshBase &boundary_mesh,
            const libMesh::DofMap &volume_dof_map,
            const libMesh::DofMap &surface_dof_map,
            unsigned int var_num,
            unsigned int comp,
            unsigned int b_var_num,
            unsigned int b_comp,
            USparseMatrix &matrix)
        {
            std::vector<libMesh::dof_id_type> dof_boundary;

            unsigned int sys_num   = volume_dof_map.sys_number();
            unsigned int b_sys_num = surface_dof_map.sys_number();

            matrix = local_sparse(volume_dof_map.n_local_dofs(), surface_dof_map.n_local_dofs(), 1);

            utopia::Write<USparseMatrix> w(matrix);
            auto rr = row_range(matrix);

            // loop through all boundary elements.
            for (const auto & b_elem : boundary_mesh.active_local_element_ptr_range())
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
                                                            b_sys_num,
                                                            b_var_num,
                                                            b_comp);
                           
                            if(rr.inside(b_dof)){
                                matrix.set(v_dof, b_dof, 1.0);
                            }
                        }
                    }
                }
            }
        }

        class Map {
        public:
            SizeType from_range_begin, from_range_end;
            SizeType to_range_begin, to_range_end;
            std::vector<SizeType> idx;

            inline SizeType from_extent() const { return from_range_end - from_range_begin; }
            inline SizeType to_extent()   const { return to_range_end   - to_range_begin; }
        };

        static void bundary_permutation_matrix_from_map(
            Map &map,
            USparseMatrix &mat)
        {
            auto n_local_dof_vol = map.to_range_end - map.to_range_begin;
            auto n_local_dofs_surf = map.idx.size();

            mat = local_sparse(n_local_dof_vol, n_local_dofs_surf, 1);

            Write<USparseMatrix> w(mat);
            for(std::size_t i = 0; i < n_local_dofs_surf; ++i) {
                mat.set(map.idx[i], i + map.from_range_begin, 1.);
            }
        }

        static void bundary_permutation_map(
            const libMesh::MeshBase &boundary_mesh,
            const libMesh::DofMap &volume_dof_map,
            const libMesh::DofMap &surface_dof_map,
            unsigned int var_num,
            unsigned int comp,
            unsigned int b_var_num,
            unsigned int b_comp,
            Map &map)
        {
            std::vector<libMesh::dof_id_type> dof_boundary;

            unsigned int sys_num   = volume_dof_map.sys_number();
            unsigned int b_sys_num = surface_dof_map.sys_number();

            map.idx.resize(surface_dof_map.n_local_dofs());
            std::fill(map.idx.begin(), map.idx.end(), 0);

            Range rr(surface_dof_map.first_dof(), surface_dof_map.last_dof() + 1);

            map.from_range_begin = rr.begin();
            map.from_range_end = rr.end();

            map.to_range_begin = volume_dof_map.first_dof();
            map.to_range_end   = volume_dof_map.last_dof() + 1;

            // loop through all boundary elements.
            for (const auto & b_elem : boundary_mesh.active_local_element_ptr_range())
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
                                                            b_sys_num,
                                                            b_var_num,
                                                            b_comp);
                           
                            if(rr.inside(b_dof)){
                                map.idx[b_dof - rr.begin()] = v_dof;
                            }
                        }
                    }
                }
            }
        }

        inline libMesh::dof_id_type n_local_elems() const
        {
            return mesh().n_active_local_elem();
        }

        inline libMesh::dof_id_type n_local_dofs() const
        {
            return n_local_dofs_;
        }

        void make_tensor_product_permutation(const int dim)
        {
            assert(permutation_);
            if(!permutation_) { return; }

            tensor_permutation_ = std::make_shared<USparseMatrix>();

            auto n_local_dofs_surf = boundary_to_volume_map.from_extent();

            assert(n_local_dofs_surf > 0);

            *tensor_permutation_ = local_sparse(
                boundary_to_volume_map.to_extent(),
                n_local_dofs_surf * dim,
                1
            );

            Write<USparseMatrix> w(*tensor_permutation_);
            for(std::size_t i = 0; i < n_local_dofs_surf; ++i) {
                for(int d = 0; d < dim; ++d) {
                    tensor_permutation_->set(
                        boundary_to_volume_map.idx[i] + d,
                        i * dim + d + boundary_to_volume_map.from_range_begin,
                        1.
                    );
                }
            }
        }

    private:
        std::shared_ptr<libMesh::MeshBase> mesh_;
        std::vector<ElementDofMap> dof_map_;
        std::vector<libMesh::dof_id_type> handle_to_element_id_;
        std::vector<FEType> fe_type_;
        std::unique_ptr<libMesh::EquationSystems> es;

        Map boundary_to_volume_map;
        std::shared_ptr<USparseMatrix> permutation_, tensor_permutation_;

        bool is_extracted_surface_;
        libMesh::dof_id_type n_local_dofs_;
    };

    // class PermutationBuilder {
    // public:
    //     //1) build permutation from surface to volume
    //     //2) build permutation from element-node dofmap to volume
    //     //3) build permutation from surface element-node dofmap to volume


       

        
    // };

    class ContactAssembler {

    };


    template<class Bound, class Collection>
    class LibMeshElemenAdapter : public moonolith::ElementAdapterBase<Bound, Collection> {
    public:
        using super = moonolith::ElementAdapterBase<Bound, Collection>;
        using super::super;

         inline const ElementDofMap &dofs() const
         {
            assert(dofs_);
            return *dofs_;
         }

        void set_dofs(const ElementDofMap * dofs)
        {
            dofs_ = dofs;
        }

    private:
        const ElementDofMap * dofs_ = nullptr;
    };  

}

namespace moonolith {

    template<class Bound>
    class ElementAdapter<Bound, utopia::LibMeshFunctionSpaceAdapter> 
    : public utopia::LibMeshElemenAdapter<Bound, utopia::LibMeshFunctionSpaceAdapter> {
    public:
        using super = utopia::LibMeshElemenAdapter<Bound, utopia::LibMeshFunctionSpaceAdapter>;
        using super::super;

    };

    template<>
    class CollectionManager<utopia::LibMeshFunctionSpaceAdapter> 
    : public utopia::LibMeshCollectionManager<utopia::LibMeshFunctionSpaceAdapter> {
    public:
        using super = utopia::LibMeshCollectionManager<utopia::LibMeshFunctionSpaceAdapter>;
        using super::super;
    };
    
}

namespace utopia {
    using LibMeshCollectionManagerT = moonolith::CollectionManager<utopia::LibMeshFunctionSpaceAdapter>;
}




#endif //UTOPIA_CONTACT_ASSEMBLER_HPP