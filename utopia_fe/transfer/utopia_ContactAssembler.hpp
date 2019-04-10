#ifndef UTOPIA_CONTACT_ASSEMBLER_HPP
#define UTOPIA_CONTACT_ASSEMBLER_HPP

#include "libmesh/mesh.h"
#include "libmesh/elem.h"
#include "libmesh/serial_mesh.h"

#include "moonolith_mesh_adapter.hpp"
#include "moonolith_polygon.hpp"
#include "moonolith_line.hpp"
#include "moonolith_check_stream.hpp"
#include "moonolith_make_unique.hpp"

#include "utopia_ElementDofMap.hpp"

#include <cassert>

namespace utopia {

    template<int Dim, class FunctionSpace>
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

        static const Elem &elem(const FunctionSpace &space, const ElementIter &e_it) 
        {
            return **e_it;
        }

        static Integer tag(const FunctionSpace &space, const ElementIter &e_it)
        {
            return (*e_it)->subdomain_id();
        }

        static Integer n_elements(const FunctionSpace &space)
        {
            return space.mesh().n_active_local_elem();
        }

        static ElementIter elements_begin(const FunctionSpace &space)
        {
            if(space.mesh().n_active_local_elem() == 0) {
                return space.mesh().local_elements_begin();
            }

            return space.mesh().active_local_elements_begin();
        }

        static ElementIter elements_end(const FunctionSpace &space)
        {
            if(space.mesh().n_active_local_elem() == 0) {
                return space.mesh().local_elements_end();
            }

            return space.mesh().active_local_elements_end();
        }

        static Integer handle(const FunctionSpace &space, const ElementIter &e_it)
        {
            return (*e_it)->id();
        }

        static bool skip(const FunctionSpace &space, const Integer &element_index)
        {
            return false;
        }

        template<class Bound>
        static void fill_bound(const FunctionSpace &space, const Integer handle, Bound &bound, const double blow_up)
        {
            const auto &mesh = space.mesh();
            const auto &e = mesh.elem(handle);
            auto dim = mesh.spatial_dimension();

            assert(dim == Dim);

            if(Dim > mesh.mesh_dimension()) {
                std::array<double, Dim> p, q;
                
                Point nn;
                compute_side_normal(Dim, e, nn);

                for(Integer i = 0; i < e.n_nodes(); ++i) {
                    
                    const auto &q = e.node_ref(i);

                    for(Integer d = 0; d < Dim; ++d) {
                        p[d] = q[d] + blow_up * nn(d);
                    }

                    bound += p.values;

                    for(Integer d = 0; d < Dim; ++d) {
                        p[d] = q[d] - blow_up * nn(d);
                    }

                    bound += p;
                }

            } else {

                std::array<double, Dim> p;
                
                for(Integer i = 0; i < n_nodes(e); ++i) {
                    const auto &q = e.node_ref(i);
                    
                    for(Integer d = 0; d < Dim; ++d) {
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

            const int dim         = space.mesh_dimension();
            const long n_elements = std::distance(begin, end);
            const auto &dof_map   = space.dof_map();
            const auto &fe_types  = space.fe_types();


            std::set<long> nodeIds;
            std::map<long, long> mapping;

            for(Iterator it = begin; it != end; ++it) {

                const libMesh::dof_id_type local_element_id = *it;
                const libMesh::dof_id_type global_element_id = space.handle_to_element_id(local_element_id);

                const libMesh::Elem *elem = space.elem(global_element_id);

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

                const libMesh::Point &p = space.node(node_id);

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

                const libMesh::Elem *elem = space.elem(global_element_id);

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

        class FEType {
        public:
            int family;
            int order;
        };

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

    private:
        std::shared_ptr<libMesh::MeshBase> mesh_;
        std::vector<ElementDofMap> dof_map_;
        std::vector<libMesh::dof_id_type> handle_to_element_id_;
        std::vector<FEType> fe_type_;
    };

    class ContactAssembler {

    };

}

#endif //UTOPIA_CONTACT_ASSEMBLER_HPP