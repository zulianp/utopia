#ifndef UTOPIA_CONTACT_ASSEMBLER_HPP
#define UTOPIA_CONTACT_ASSEMBLER_HPP

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

            bundary_permutation_matrix(
            *b_mesh,
            dof_map,
            sys.get_dof_map(),
            var_num,
            0,
            0,
            0,
            *permutation_);

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

        inline libMesh::dof_id_type n_local_dofs()
        {
            return n_local_dofs_;
        }

    private:
        std::shared_ptr<libMesh::MeshBase> mesh_;
        std::vector<ElementDofMap> dof_map_;
        std::vector<libMesh::dof_id_type> handle_to_element_id_;
        std::vector<FEType> fe_type_;
        std::unique_ptr<libMesh::EquationSystems> es;

        std::shared_ptr<USparseMatrix> permutation_;

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

namespace utopia {

        template<typename T, int Dim>
        inline void make(const libMesh::Point &p, moonolith::Vector<T, Dim> &q)
        {
            for(int i = 0; i < Dim; ++i) {
                q[i] = p(i);
            }
        }

        template<typename T, int Dim>
        inline void make(const libMesh::Elem &elem, moonolith::Line<T, Dim> &poly)
        {   
            make(elem.node_ref(0), poly.p0);
            make(elem.node_ref(1), poly.p1);
        }


        template<typename T, int Dim>
        inline void make_non_affine(const libMesh::Elem &elem, moonolith::Storage<moonolith::Vector<T, Dim>> &poly_line)
        {   
            const int n_nodes = elem.n_nodes();

            poly_line.resize(n_nodes);

            if(n_nodes == 2) {
                //P1
                make(elem.node_ref(0), poly_line[0]);
                make(elem.node_ref(1), poly_line[1]);
            } else if(n_nodes == 3) {
                //P2
                make(elem.node_ref(0), poly_line[0]);
                make(elem.node_ref(2), poly_line[1]);
                make(elem.node_ref(1), poly_line[2]);
            } else if(n_nodes == 4) {
                //P3
                make(elem.node_ref(0), poly_line[0]);
                make(elem.node_ref(2), poly_line[1]);
                make(elem.node_ref(3), poly_line[2]);
                make(elem.node_ref(1), poly_line[3]);
            } else {
                assert(false);
            }

        }

        template<typename T, int Dim>
        inline void make_non_affine(const libMesh::Elem &elem, moonolith::PolyLine<T, Dim> &poly_line)
        {
            make_non_affine(elem, poly_line.points);
        }

        template<typename T, int Dim>
        inline void make(const libMesh::Elem &elem, moonolith::Polygon<T, Dim> &poly)
        {
           assert(is_tri(elem.type()) || is_quad(elem.type()));
           auto n_nodes = is_tri(elem.type()) ? 3 : 4;

           poly.resize(n_nodes);

           for(auto i = 0; i < n_nodes; ++i) {
               const auto &p = elem.node_ref(i);
               auto &q = poly.points[i];
               make(p, q);
           }
        }


        template<typename T, int Dim>
        inline void make_non_affine(const libMesh::Elem &elem, moonolith::Polygon<T, Dim> &poly)
        {
           assert(is_tri(elem.type()) || is_quad(elem.type()));
           auto n_corner_nodes = is_tri(elem.type()) ? 3 : 4;
           auto n_nodes = elem.n_nodes();

           if(n_nodes == n_corner_nodes) {
                //fallback to affine case
                return make(elem, poly);
           }

           //P2 only
           assert(n_nodes = n_corner_nodes * 2);

           poly.resize(n_nodes);

           for(auto i = 0; i < n_nodes; i+=2) {
               const auto &p = elem.node_ref(i);
               const auto &q = elem.node_ref(i + n_corner_nodes);
               make(p, poly.points[i]);
               make(q, poly.points[i + 1]);
           }
        }

        template<typename T>
        inline void make_tetrahedron(const libMesh::Elem &elem, moonolith::Polyhedron<T> &poly)
        {
            static const std::size_t n_nodes = 4;

            assert(elem.dim() == 3);
            assert(is_tet(elem.type()));

            poly.el_ptr.resize(n_nodes + 1);
            poly.el_index.resize(12);
            poly.points.resize(n_nodes);

            poly.el_ptr[0] = 0;
            poly.el_ptr[1] = 3;
            poly.el_ptr[2] = 6;
            poly.el_ptr[3] = 9;
            poly.el_ptr[4] = 12;

            for(std::size_t i = 0; i < n_nodes; ++i) {
                const  auto &p = elem.node_ref(i);
                poly.points[i].x = p(0);
                poly.points[i].y = p(1);
                poly.points[i].z = p(2);
            }

            poly.el_index[0] = 0;
            poly.el_index[1] = 2;
            poly.el_index[2] = 1;
       
            poly.el_index[3] = 0;
            poly.el_index[4] = 3;
            poly.el_index[5] = 2;
       
            poly.el_index[6] = 0;
            poly.el_index[7] = 1;
            poly.el_index[8] = 3;
       
            poly.el_index[9] = 1;
            poly.el_index[10] = 2;
            poly.el_index[11] = 3;   

            poly.type = moonolith::Polyhedron<T>::TET;     
        }

        template<typename T>
        void make_pyramid(const libMesh::Elem &elem, moonolith::Polyhedron<T> &poly)
        {
            static const std::size_t n_nodes = 5;
            assert(elem.dim() == 3);
            assert(is_pyramid(elem.type()));

            poly.el_ptr.resize(5+1);
            poly.el_index.resize(16);
            poly.points.resize(n_nodes);

            for(std::size_t i = 0; i < n_nodes; ++i) {
                const  auto &p = elem.node_ref(i);
                poly.points[i].x = p(0);
                poly.points[i].y = p(1);
                poly.points[i].z = p(2);
            }

            poly.el_ptr[0] = 0;
            poly.el_ptr[1] = 3;
            poly.el_ptr[2] = 6;
            poly.el_ptr[3] = 9;
            poly.el_ptr[4] = 12;
            poly.el_ptr[5] = 16;

            //face 0
            poly.el_index[0] = 0;
            poly.el_index[1] = 1;
            poly.el_index[2] = 4;

            //face 1
            poly.el_index[3] = 1;
            poly.el_index[4] = 2;
            poly.el_index[5] = 4;

            //face 2
            poly.el_index[6] = 3;
            poly.el_index[7] = 0;
            poly.el_index[8] = 4;

            //face 3
            poly.el_index[9]  = 2;
            poly.el_index[10] = 3;
            poly.el_index[11] = 4;

            //face 4
            poly.el_index[12] = 0;
            poly.el_index[13] = 3;
            poly.el_index[14] = 2;
            poly.el_index[15] = 1;

            poly.type = moonolith::Polyhedron<T>::UNSTRUCTURED;
        }

        template<typename T>
        void make_prism(const libMesh::Elem &elem, moonolith::Polyhedron<T> &poly)
        {
            assert(elem.dim() == 3);
            assert(is_prism(elem.type()));

            static const std::size_t n_nodes = 6;

            poly.el_ptr.resize(5+1);
            poly.el_index.resize(18);
            poly.points.resize(n_nodes);

            for(std::size_t i = 0; i < n_nodes; ++i) {
                const  auto &p = elem.node_ref(i);
                poly.points[i].x = p(0);
                poly.points[i].y = p(1);
                poly.points[i].z = p(2);
            }

            poly.el_ptr[0] = 0;
            poly.el_ptr[1] = 3;
            poly.el_ptr[2] = 6;
            poly.el_ptr[3] = 10;
            poly.el_ptr[4] = 14;
            poly.el_ptr[5] = 18;

            //face 0
            poly.el_index[0] = 0;
            poly.el_index[1] = 1;
            poly.el_index[2] = 2;

            //face 1
            poly.el_index[3] = 4;
            poly.el_index[4] = 3;
            poly.el_index[5] = 5;

            //face 2
            poly.el_index[6] = 0;
            poly.el_index[7] = 2;
            poly.el_index[8] = 5;
            poly.el_index[9] = 3;

            //face 3
            poly.el_index[10] = 1;
            poly.el_index[11] = 4;
            poly.el_index[12] = 5;
            poly.el_index[13] = 2;

            //face 4
            poly.el_index[14] = 0;
            poly.el_index[15] = 1;
            poly.el_index[16] = 4;
            poly.el_index[17] = 3;

            poly.type = moonolith::Polyhedron<T>::UNSTRUCTURED;
        }

        template<typename T>
        void make_hex(const libMesh::Elem &elem, moonolith::Polyhedron<T> &poly)
        {
            assert(elem.dim() == 3);
            assert(is_hex(elem.type()));

            static const std::size_t n_nodes = 8;
            poly.el_ptr.resize(6+1);
            poly.el_index.resize(6*4);
            poly.points.resize(n_nodes);

            for(std::size_t i = 0; i < n_nodes; ++i) {
                const  auto &p = elem.node_ref(i);
                poly.points[i].x = p(0);
                poly.points[i].y = p(1);
                poly.points[i].z = p(2);
            }

            poly.el_ptr[0] = 0;
            poly.el_ptr[1] = 4;
            poly.el_ptr[2] = 8;
            poly.el_ptr[3] = 12;
            poly.el_ptr[4] = 16;
            poly.el_ptr[5] = 20;
            poly.el_ptr[6] = 24;

            //face 0
            poly.el_index[0] = 0;
            poly.el_index[1] = 1;
            poly.el_index[2] = 5;
            poly.el_index[3] = 4;

            //face 1
            poly.el_index[4] = 1;
            poly.el_index[5] = 2;
            poly.el_index[6] = 6;
            poly.el_index[7] = 5;

            //face 2
            poly.el_index[8]  = 3;
            poly.el_index[9]  = 7;
            poly.el_index[10] = 6;
            poly.el_index[11] = 2;

            //face 3
            poly.el_index[12] = 0;
            poly.el_index[13] = 4;
            poly.el_index[14] = 7;
            poly.el_index[15] = 3;

            //face 4
            poly.el_index[16] = 2;
            poly.el_index[17] = 1;
            poly.el_index[18] = 0;
            poly.el_index[19] = 3;

            //face 5
            poly.el_index[20] = 6;
            poly.el_index[21] = 7;
            poly.el_index[22] = 4;
            poly.el_index[23] = 5;

            poly.type = moonolith::Polyhedron<T>::HEX;
        }

        template<typename T>
        inline void make(const libMesh::Elem &elem, moonolith::Polyhedron<T> &poly)
        {
            if(is_tet(elem.type())) {
                make_tetrahedron(elem, poly);
                return;
            }

            if(is_hex(elem.type())) {
                make_hex(elem, poly);
                return;
            }

            if(is_pyramid(elem.type())) {
                make_pyramid(elem, poly);
                return;
            }

            if(is_prism(elem.type())) {
                make_prism(elem, poly);
                return;
            }

            assert(false);
        }


    template<typename T, int Dim>
    inline void normal(const libMesh::Elem &elem, moonolith::Vector<T, Dim> &nn)
    {
        using namespace libMesh;
        Point o, u, v, n;

        if(Dim == 2) {
            assert(elem.n_nodes() >= 2);
            o = elem.point(0);
            u = elem.point(1);
            u -= o;
            n(0) =  u(1);
            n(1) = -u(0);

        } else {
            assert(Dim == 3);
            o = elem.point(0);
            u = elem.point(1);
            v = elem.point(2);
            u -= o;
            v -= o;
            n = u.cross(v);
        }

        n *= 1./n.norm();

        make(n, nn);
    }

    template<typename T, int Dim>
    inline void convert(
        const libMesh::QGauss &q_in,
        const moonolith::Vector<T, Dim> &point_shift,
        const T &point_rescale,
        const T &weight_rescale,
        moonolith::Quadrature<T, Dim> &q_out)
    {
        const auto &p_in = q_in.get_points();
        const auto &w_in = q_in.get_weights();

        const std::size_t n_qp = q_in.n_points();
        
        q_out.resize(n_qp);

        for(std::size_t k = 0; k < n_qp; ++k) {
            q_out.weights[k] = w_in[k] * weight_rescale;

            const auto &pk_in = p_in[k];
            auto &pk_out = q_out.points[k];

            for(int i = 0; i < Dim; ++i) {
                pk_out[i] = (pk_in(i) + point_shift[i]) * point_rescale;
            }
        }
    }

    template<typename T, int Dim>
    inline void convert(
        const moonolith::Quadrature<T, Dim> &q_in,
        const moonolith::Vector<T, Dim> &point_shift,
        const T &point_rescale,
        const T &weight_rescale,
        QMortar &q_out)
    {
        const auto &p_in = q_in.points;
        const auto &w_in = q_in.weights;

        auto &p_out = q_out.get_points();
        auto &w_out = q_out.get_weights();

        const std::size_t n_qp = q_in.n_points();
        
        q_out.resize(n_qp);

        for(std::size_t k = 0; k < n_qp; ++k) {
            w_out[k] = w_in[k] * weight_rescale;

            const auto &pk_in = p_in[k];
            auto &pk_out = p_out[k];

            for(int i = 0; i < Dim; ++i) {
                pk_out(i) = (pk_in[i] + point_shift[i]) * point_rescale;
            }
        }
    }

    template<typename T, int Dim>
    inline void convert(
        const moonolith::Quadrature<T, Dim> &q_in,
        const T &weight_rescale,
        QMortar &q_out)
    {
        static const moonolith::Vector<T, Dim> zero;
        return convert(q_in, zero, 1., weight_rescale, q_out);
    }

    // inline void make_line_transform(
    //     const libMesh::Elem &elem,
    //     moonolith::AffineTransform<double, 1, 2> &trafo)
    // {
    //     const auto &q0 = elem.node_ref(0);
    //     const auto &q1 = elem.node_ref(1);

    //     moonolith::Vector<double, 2> p0, p1;

    //     p0.x = q1(0);
    //     p0.y = q0(1);

    //     p1.x = q1(0);
    //     p1.y = q1(1);

    //     make(p0, p1, trafo);
    // }

    // inline void make_triangle_transform(
    //     const libMesh::Elem &elem,
    //     moonolith::AffineTransform<double, 2, 3> &trafo)
    // {
    //     const auto &q0 = elem.node_ref(0);
    //     const auto &q1 = elem.node_ref(1);
    //     const auto &q2 = elem.node_ref(2);

    //     moonolith::Vector<double, 3> p0, p1, p2;

    //     p0.x = q0(0);
    //     p0.y = q0(1);
    //     p0.z = q0(2);

    //     p1.x = q1(0);
    //     p1.y = q1(1);
    //     p1.z = q1(2);

    //     p2.x = q2(0);
    //     p2.y = q2(1);
    //     p2.z = q2(2);

    //     make(p0, p1, p2, trafo);
    // }

    inline void make_transform(
        const libMesh::Elem &elem,
        std::shared_ptr<moonolith::Transform<double, 2, 3>> &trafo)
    {
        trafo = std::make_shared<Transform2>(elem);
    }   

    inline void make_transform(
        const libMesh::Elem &elem,
        std::shared_ptr<moonolith::Transform<double, 1, 2>> &trafo)
    {
        trafo = std::make_shared<Transform1>(elem);
    }   

    inline void make_transform(const libMesh::Elem &elem, moonolith::AffineTransform<double, 2, 2> &trafo)
    {
        libMesh::Point p0(0.0, 0.0, 0.0);
        libMesh::Point p1(1.0, 0.0, 0.0);
        libMesh::Point p2(0.0, 1.0, 0.0);

        // figure out the ordering of the map
        p0 = libMesh::FE<2, libMesh::LAGRANGE>::map(&elem, p0);
        p1 = libMesh::FE<2, libMesh::LAGRANGE>::map(&elem, p1);
        p2 = libMesh::FE<2, libMesh::LAGRANGE>::map(&elem, p2);

        moonolith::Vector<double, 2> q0, q1, q2;

        q0.x = p0(0);
        q0.y = p0(1);

        q1.x = p1(0);
        q1.y = p1(1);

        q2.x = p2(0);
        q2.y = p2(1);

        moonolith::make(q0, q1, q2, trafo);
    }

    inline void make_transform(const libMesh::Elem &elem, moonolith::AffineTransform<double, 3, 3> &trafo)
    {
        libMesh::Point p0(0.0, 0.0, 0.0);
        libMesh::Point p1(1.0, 0.0, 0.0);
        libMesh::Point p2(0.0, 1.0, 0.0);
        libMesh::Point p3(0.0, 0.0, 1.0);

        // figure out the ordering of the map
        p0 = libMesh::FE<3, libMesh::LAGRANGE>::map(&elem, p0);
        p1 = libMesh::FE<3, libMesh::LAGRANGE>::map(&elem, p1);
        p2 = libMesh::FE<3, libMesh::LAGRANGE>::map(&elem, p2);
        p3 = libMesh::FE<3, libMesh::LAGRANGE>::map(&elem, p3);

        moonolith::Vector<double, 3> q0, q1, q2, q3;

        q0.x = p0(0);
        q0.y = p0(1);
        q0.z = p0(2);

        q1.x = p1(0);
        q1.y = p1(1);
        q1.z = p1(2);

        q2.x = p2(0);
        q2.y = p2(1);
        q2.z = p2(2);

        q3.x = p3(0);
        q3.y = p3(1);
        q3.z = p3(2);

        moonolith::make(q0, q1, q2, q3, trafo);
    }

    inline void make_transform(const libMesh::Elem &elem, moonolith::AffineTransform<double, 2, 3> &trafo)
    {
        libMesh::Point p0(0.0, 0.0);
        libMesh::Point p1(1.0, 0.0);
        libMesh::Point p2(0.0, 1.0);

        // figure out the ordering of the map
        p0 = libMesh::FE<2, libMesh::LAGRANGE>::map(&elem, p0);
        p1 = libMesh::FE<2, libMesh::LAGRANGE>::map(&elem, p1);
        p2 = libMesh::FE<2, libMesh::LAGRANGE>::map(&elem, p2);

        moonolith::Vector<double, 3> q0, q1, q2;

        q0.x = p0(0);
        q0.y = p0(1);
        q0.z = p0(2);

        q1.x = p1(0);
        q1.y = p1(1);
        q1.z = p1(2);

        q2.x = p2(0);
        q2.y = p2(1);
        q2.z = p2(2);

        moonolith::make(q0, q1, q2, trafo);
    }

    inline void make_transform(const libMesh::Elem &elem, moonolith::AffineTransform<double, 1, 2> &trafo)
    {
        libMesh::Point p0(-1.0, 0.0);
        libMesh::Point p1(1.0, 0.0);

        // figure out the ordering of the map
        p0 = libMesh::FE<1, libMesh::LAGRANGE>::map(&elem, p0);
        p1 = libMesh::FE<1, libMesh::LAGRANGE>::map(&elem, p1);

        moonolith::Vector<double, 2> q0, q1;

        q0.x = p0(0);
        q0.y = p0(1);

        q1.x = p1(0);
        q1.y = p1(1);

        moonolith::make(q0, q1, trafo);
    }


}



#endif //UTOPIA_CONTACT_ASSEMBLER_HPP