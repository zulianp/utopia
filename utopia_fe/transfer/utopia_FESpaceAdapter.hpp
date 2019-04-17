#ifndef UTOPIA_FE_SPACE_ADAPTER_HPP
#define UTOPIA_FE_SPACE_ADAPTER_HPP

#include "utopia_ElementDofMap.hpp"
#include "moonolith_communicator.hpp"

#include "libmesh/serial_mesh.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"

#include <memory>
#include <vector>
#include <utility>

namespace libMesh {
    class MeshBase;
    class DofMap;
}

namespace utopia {
    class FESpaceAdapter {
    public:
        inline explicit FESpaceAdapter(const moonolith::Communicator &comm) : comm(comm) {}

        FESpaceAdapter(
            const std::shared_ptr<libMesh::MeshBase> &master_slave,
            const std::shared_ptr<libMesh::DofMap>   &original_dofmap,
            const unsigned int var_num,  const std::vector< std::pair<int, int>> &tags);

        inline std::shared_ptr<libMesh::MeshBase> &mesh()
        {
            return mesh_;
        }

        inline const std::shared_ptr<libMesh::MeshBase> &mesh() const
        {
            return mesh_;
        }

        inline long n_elements() const
        {
            long ret=0;
            ret += mesh_->n_elem();
            return ret;
        }

        inline std::vector<ElementDofMap> &dof_map()
        {

            return dof_maps_;
        }

        inline const std::vector<ElementDofMap> &dof_map() const
        {

            return dof_maps_;
        }

        inline std::vector<long> &variable_number()
        {

            return var_number_;
        }

        inline const std::vector<long> &variable_number() const
        {
            return var_number_;
        }

        inline std::vector<long> &variable_order()
        {

            return var_order_;
        }

        inline const std::vector<long> &variable_order() const
        {

            return var_order_;
        }

        inline std::vector<long> &variable_type()
        {

            return var_type_;
        }

        inline const std::vector<long> &variable_type() const
        {
            return var_type_;
        }

        inline std::vector<ElementDofMap> &subdomain_id()
        {

            return subdomain_id_;
        }

        inline const std::vector<ElementDofMap> &subdomain_id() const
        {
            return subdomain_id_;
        }

        inline std::vector<ElementDofMap> &side_set_id()
        {

            return side_set_id_;
        }

        inline const std::vector<ElementDofMap> & side_set_id() const
        {
            return side_set_id_;
        }

        inline std::vector<ElementDofMap> &face_set_id_global()
        {

            return face_set_id_global_;
        }

        inline const std::vector<ElementDofMap> & face_set_id_global() const
        {
            return face_set_id_global_;
        }

        inline std::vector<ElementDofMap> & n_face_nodes()
        {

            return n_face_nodes_;
        }

        inline const std::vector<ElementDofMap> & n_face_nodes() const
        {
            return n_face_nodes_;
        }

        inline std::vector<ElementDofMap> &side_set_id_tag()
        {

            return side_set_id_tag_;
        }

        inline const std::vector<ElementDofMap> &side_set_id_tag() const
        {
            return side_set_id_tag_;
        }


        inline std::vector<moonolith::Integer> &ownershipRangesFaceID()
        {

            return ownershipRangesFaceID_;
        }

        inline const std::vector<moonolith::Integer> &ownershipRangesFaceID() const
        {
            return ownershipRangesFaceID_;
        }

        inline std::vector<libMesh::dof_id_type> &handle_to_element_id()
        {
            return handle_to_element_id_;
        }

        inline const std::vector<libMesh::dof_id_type> &handle_to_element_id() const
        {
            return handle_to_element_id_;
        }

    private:
        moonolith::Communicator comm;
        std::shared_ptr< libMesh::MeshBase> mesh_;
        std::vector<ElementDofMap> dof_maps_;
        std::vector<libMesh::dof_id_type> handle_to_element_id_;

        std::vector<long> var_number_;
        std::vector<long> var_type_;
        std::vector<long> var_order_;




        std::vector<ElementDofMap> side_set_id_;
        std::vector<ElementDofMap> face_set_id_global_;
           std::vector<ElementDofMap> subdomain_id_;
        std::vector<ElementDofMap> side_set_id_tag_;
        std::vector<ElementDofMap> n_face_nodes_;
        std::vector<moonolith::Integer> ownershipRangesFaceID_;

        static bool is_tagged_contact_boundary(
            const libMesh::MeshBase &mesh,
            const libMesh::Elem *elem,
            const int side_elem,
            const std::vector< std::pair<int, int> > &tags);


        static void copy_global_dofs(
            libMesh::MeshBase &mesh,
            const std::shared_ptr<libMesh::DofMap>  &original_dofmap,
            const unsigned int  var_num,
            std::vector<ElementDofMap> &dof_map,
            std::vector<long> &variable_type,
            std::vector<long> &variable_number,
            std::vector<ElementDofMap> &subdomain_id,
            std::vector<ElementDofMap> &side_set_id,
            std::vector<ElementDofMap> &side_set_id_tag,
            std::vector<ElementDofMap> &face_set_id_global,
            std::vector<moonolith::Integer> &ownershipRangesFaceID,
            std::vector<libMesh::dof_id_type> &handle_to_element_id,
            const std::vector< std::pair<int, int> > &tags);


        static void copy_var_order(libMesh::DofMap &dofmap, std::vector<long> &variable_order);
    };

    template<class Iterator>
    static void write_space(
                            const Iterator &begin,
                            const Iterator &end,
                            libMesh::MeshBase &space,
                            const std::vector<ElementDofMap> &dof_map,
                            const std::vector<long> &variable_number,
                            const std::vector<long> &variable_order,
                            const std::vector<ElementDofMap> &subdomain_id,
                            const std::vector<ElementDofMap> &side_set_id,
                            const std::vector<ElementDofMap> &face_set_id_global,
                            const std::vector<libMesh::dof_id_type> &handle_to_element_id,
                            moonolith::OutputStream &os)
    {
        const int dim 		  = space.mesh_dimension();
        const long n_elements = std::distance(begin, end);

        std::set<long> nodeIds;
        std::map<long, long> mapping;
        std::vector<libMesh::dof_id_type> dof_array;

        for(Iterator it = begin; it != end; ++it) {
            //ID_FIX
            const libMesh::dof_id_type global_element_id = handle_to_element_id[*it];
            // const libMesh::Elem *elem = space.elem(*it);
            const libMesh::Elem *elem = space.elem(global_element_id);

            for(libMesh::dof_id_type j = 0; j != elem->n_nodes(); ++j) {

                nodeIds.insert(elem->node(j));


            }
        }

        long n_nodes = nodeIds.size();

        // Estimate for allocation
        os.request_space( (n_elements * 8 + n_nodes * dim) * (sizeof(double) + sizeof(long)) );

        //WRITE 1
        os << dim;

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

            for(int i = 0; i < dim; ++i) {

                //WRITE 3
                os << p(i);

            }
        }

        std::vector<libMesh::dof_id_type> indices_vector;



        for(Iterator it = begin; it != end; ++it) {

            //			const int k = *it;
            //ID_FIX
            const libMesh::dof_id_type local_element_id = *it;
            const libMesh::dof_id_type global_element_id = handle_to_element_id[*it];
            const libMesh::Elem *elem = space.elem(global_element_id);
            // const libMesh::Elem *elem = space.elem(*it);

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

            //			bool  size=true;

            int volume_tag;

            volume_tag=subdomain_id[local_element_id].global.at(0);

            os << volume_tag;

            int side_set_tag;

            //			int face_id;

            //			bool check_side_id_one = true;
            side_set_tag=side_set_id[local_element_id].global.at(0);

            os << side_set_tag;
            os << face_set_id_global.at(local_element_id);
        }
        //
        //

        //WRITE 11
        os << variable_number.at(0);

        //WRITE 13
        os << variable_order.at(0);


    }

    template<class Iterator>
    static void write_element_selection(
                                        const Iterator &begin,
                                        const Iterator &end,
                                        const FESpaceAdapter &fespace,
                                        moonolith::OutputStream &os)
    {
        write_space(
            begin,
            end,
            *fespace.mesh(),
            fespace.dof_map(),
            fespace.variable_number(),
            fespace.variable_number(),
            fespace.subdomain_id(),
            fespace.side_set_id(),
            fespace.face_set_id_global(),
            fespace.handle_to_element_id(),
            os);
    }


    static void read_space(moonolith::InputStream &is,
                           std::shared_ptr<libMesh::MeshBase> & space,
                           std::vector<ElementDofMap> &dof_map,
                           std::vector<long> &variable_number,
                           std::vector<long> &variable_order,
                           std::vector<ElementDofMap> &subdomain_id,
                           std::vector<ElementDofMap> &side_set_id,
                           std::vector<ElementDofMap> &face_set_id_global,
                           std::vector<libMesh::dof_id_type> &handle_to_element_id,
                           const libMesh::Parallel::Communicator &comm)
    {
        using namespace std;

        //READ 1
        int dim;
        is >> dim;

        //READ 2
        long n_nodes;
        is >> n_nodes;

        //READ 6
        long n_elements;
        is >> n_elements;

        auto mesh_ptr = std::make_shared<libMesh::SerialMesh>(comm, dim);

        mesh_ptr->reserve_nodes(n_nodes);

        for (long iii = 0; iii != n_nodes; ++iii) {

            libMesh::Point p;

            for(int j = 0; j < dim; ++j) {
                //READ 3
                is >> p(j);
            }

            mesh_ptr->add_point(p);
        }



        dof_map.resize(n_elements);

        subdomain_id.resize(n_elements);

        side_set_id.resize(n_elements);

        face_set_id_global.resize(n_elements);

        face_set_id_global.resize(n_elements);

        handle_to_element_id.resize(n_elements);


        for(long i = 0; i !=n_elements; ++i) {

            //READ 7

            int type, e_n_nodes;

            is >> type >> e_n_nodes;

            //std::cout<<"e_n_nodes_read = "<<e_n_nodes<<std::endl;

            auto elem =  libMesh::Elem::build(libMesh::ElemType(type)).release();
            handle_to_element_id[i] = i;

            //std::cout<<"n_side_read ="<< elem->n_sides()<<std::endl;


            int index;

            for (int ii = 0; ii != e_n_nodes; ++ii) {

                //READ 8
                is >> index;
                //std::cout<<"index = "<<index<<std::endl;
                elem->set_node(ii) = & mesh_ptr->node(index);

            }


            //READ 9
            is >> dof_map.at(i);
            //std::cout<< "dof_map_read = "<<dof_map[i].global.at(0)<<std::endl;

            int volume_tag, side_set_tag;
            //			int face_id;

            //			bool on_boundary=false;
            //std::cout<<"read n_elements = "<<n_elements<<std::endl;


            is >> volume_tag;

            //std::cout<<" read volume role = "<< volume_tag <<std::endl;

            subdomain_id[i].global.insert(subdomain_id[i].global.end(),volume_tag);

            is >> side_set_tag;

            is >> face_set_id_global.at(i);

            //std::cout <<"read value"<< face_set_id_global[i].global.at(0)<<std::endl;

            side_set_id[i].global.insert(side_set_id[i].global.end(),side_set_tag);

            mesh_ptr->add_elem(elem);

            libmesh_assert(elem);

        }

        //READ 11
        variable_number.resize(1);
        is >> variable_number.at(0);

        //READ 12
        variable_order.resize(1);
        is >> variable_order.at(0);


        //!!!! dummy parameters
        space = mesh_ptr;

    }

    static void read_spaces(moonolith::InputStream &is, FESpaceAdapter &utopiamesh, const libMesh::Parallel::Communicator &comm_mesh)
    {
        read_space(
            is,
            utopiamesh.mesh(),
            utopiamesh.dof_map(),
            utopiamesh.variable_number(),
            utopiamesh.variable_order(),
            utopiamesh.subdomain_id(),
            utopiamesh.side_set_id(),
            utopiamesh.face_set_id_global(),
            utopiamesh.handle_to_element_id(),
            comm_mesh);
    }

}


#endif //UTOPIA_FE_SPACE_ADAPTER_HPP
