#ifndef UTOPIA_FE_SPACES_R_ADAPTER_HPP
#define UTOPIA_FE_SPACES_R_ADAPTER_HPP 


#include "utopia_copy_dofmap.hpp"
#include "utopia_ElementDofMap.hpp"

#include "moonolith_communicator.hpp"
#include "moonolith_check_stream.hpp"

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

    class FESpacesRAdapter {
    public:

      inline FESpacesRAdapter(const moonolith::Communicator &comm) : comm(comm){}

      FESpacesRAdapter(
       const std::shared_ptr<libMesh::MeshBase> &master,
       const std::shared_ptr<libMesh::MeshBase> &slave,
       const std::shared_ptr<libMesh::DofMap> &dof_map_master,
       const std::shared_ptr<libMesh::DofMap> &dof_map_slave,
       const std::shared_ptr<libMesh::DofMap> &dof_map_reverse_master,
       const std::shared_ptr<libMesh::DofMap> &dof_map_reverse_slave,
       const unsigned int &from_var_num,
       const unsigned int &to_var_num,
       const unsigned int &from_var_num_r,
       const unsigned int &to_var_num_r);


      inline std::vector< std::shared_ptr<libMesh::MeshBase> > &spaces()
      {
            return spaces_;
        }

    inline const std::vector< std::shared_ptr<libMesh::MeshBase> > &spaces() const
    {
        return spaces_;
    }

    inline long n_elements() const
    {
        long ret = 0;
        for(auto s : spaces_) {
            if(s) {
                ret += s->n_elem();
            }
        }

        return ret;
    }

    inline std::vector<ElementDofMap> &dof_map(const int i)
    {
        assert(i < 2);
        assert(i >= 0);
        return dof_maps_[i];
    }

    inline const std::vector<ElementDofMap> &dof_map(const int i) const
    {
        assert(i < 2);
        assert(i >= 0);
        return dof_maps_[i];
    }

    inline std::vector<ElementDofMap> &dof_map_reverse(const int i)
    {
        assert(i < 2);
        assert(i >= 0);
        return dof_maps_reverse_[i];
    }

    inline const std::vector<ElementDofMap> &dof_map_reverse(const int i) const
    {
        assert(i < 2);
        assert(i >= 0);
        return dof_maps_reverse_[i];
    }

    inline std::vector<ElementDofMap> &variable_number(const int i)
    {
        assert(i < 2);
        assert(i >= 0);
        return var_number_[i];
    }

    inline const std::vector<ElementDofMap> &variable_number(const int i) const
    {
        assert(i < 2);
        assert(i >= 0);
        return var_number_[i];
    }

    inline std::vector<ElementDofMap> &variable_order(const int i)
    {
        assert(i < 2);
        assert(i >= 0);
        return var_order_[i];
    }

    inline const std::vector<ElementDofMap> &variable_order(const int i) const
    {
        assert(i < 2);
        assert(i >= 0);
        return var_order_[i];
    }


    inline std::vector<ElementDofMap> &variable_type(const int i)
    {
        assert(i < 2);
        assert(i >= 0);
        return var_type_[i];
    }



    inline const std::vector<ElementDofMap> &variable_type(const int i) const
    {
        assert(i < 2);
        assert(i >= 0);
        return var_type_[i];
    }

    inline std::vector<libMesh::dof_id_type> &handle_to_element_id(const int i)
    {
        assert(i < 2);
        assert(i >= 0);
        return handle_to_element_id_[i];
    }

    inline const std::vector<libMesh::dof_id_type> &handle_to_element_id(const int i) const
    {
        assert(i < 2);
        assert(i >= 0);
        return handle_to_element_id_[i];
    }


private:
    moonolith::Communicator comm;
    std::vector<std::shared_ptr<libMesh::MeshBase> > spaces_;
    std::vector<ElementDofMap> dof_maps_[2];
    std::vector<ElementDofMap> dof_maps_reverse_[2];
    std::vector<ElementDofMap> var_number_[2];
    std::vector<ElementDofMap> var_order_[2];
    std::vector<ElementDofMap> var_type_[2];
    std::vector<libMesh::dof_id_type> handle_to_element_id_[2];        
};

template<class Iterator>
static void write_space(
    const Iterator &begin, 
    const Iterator &end,
    libMesh::MeshBase &space,
    const std::vector<ElementDofMap> &dof_map,
    const std::vector<ElementDofMap> &dof_map_reverse,
    const std::vector<ElementDofMap> &variable_order,
    const std::vector<libMesh::dof_id_type> &handle_to_element_id,
    const int role, 
    moonolith::OutputStream &os)
{
    const int dim         = space.mesh_dimension();
    const long n_elements = std::distance(begin, end);

    std::set<long> nodeIds;
    std::map<long, long> mapping;

    for(Iterator it = begin; it != end; ++it) {
        const libMesh::dof_id_type local_element_id = *it;
        const libMesh::dof_id_type global_element_id = handle_to_element_id[local_element_id];
        const libMesh::Elem *elem = space.elem(global_element_id);

        for(libMesh::dof_id_type j = 0; j != elem->n_nodes(); ++j) {

            nodeIds.insert(elem->node(j));
        }
    }

    long n_nodes = nodeIds.size();

    // Estimate for allocation
    os.request_space( (n_elements * 8 + n_nodes * dim) * (sizeof(double) + sizeof(long)) );

    //WRITE 1
    os << dim << role;


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

    CHECK_STREAM_WRITE_BEGIN("elements", os);

    for(Iterator it = begin; it != end; ++it) {
     const libMesh::dof_id_type local_element_id = *it;
     const libMesh::dof_id_type global_element_id = handle_to_element_id[local_element_id];

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
    assert(!dof_map_reverse.at(local_element_id).empty());

    os << dof_map.at(local_element_id);
    os << dof_map_reverse.at(local_element_id);
}

CHECK_STREAM_WRITE_END("elements", os);

    //WRITE 10
    //        os << variable_number.at(0);

    //WRITE 11
os << variable_order.at(0);



}

template<class Iterator>
static void write_element_selection(
    const Iterator &begin, 
    const Iterator &end, 
    const FESpacesRAdapter &spaces, 
    moonolith::OutputStream &os)
{
    if(spaces.spaces().empty()){
        assert(false);
        return;
    }

    auto m = spaces.spaces()[0];
    std::shared_ptr<libMesh::MeshBase> s = nullptr;

    if(spaces.spaces().size()>1) {
        s=spaces.spaces()[1];
    }

    std::vector<long> master_selection;
    std::vector<long> slave_selection;

    bool met_slave_selection = false;

    const libMesh::dof_id_type n_elem = spaces.handle_to_element_id(0).size();

    for(Iterator it = begin; it != end; ++it) {
        int index = *it;

        if(m && index >= n_elem) {
            index -= n_elem;
            slave_selection.push_back(index);
        } else if(!m) {
            met_slave_selection = true;
            slave_selection.push_back(index);
        } else {
            assert(!met_slave_selection);
            assert(index < n_elem);
            master_selection.push_back(index);
        }
    }

    const bool has_master = !master_selection.empty();
    const bool has_slave  = !slave_selection.empty();

    os << has_master << has_slave;

    if(has_master) {
        write_space(
            master_selection.begin(), 
            master_selection.end(), 
            *m, spaces.dof_map(0), 
            spaces.dof_map_reverse(0),
            spaces.variable_order(0), 
            spaces.handle_to_element_id(0),
            0, 
            os);
    }

    if(has_slave) {
        write_space(
            slave_selection.begin(),
            slave_selection.end(),
            *s, spaces.dof_map(1),
            spaces.dof_map_reverse(1),
            spaces.variable_order(1),
            spaces.handle_to_element_id(1),
            1,
            os);
    }
}

static void read_space(
    moonolith::InputStream &is, 
    std::shared_ptr<libMesh::MeshBase> & space,
    std::vector<ElementDofMap> &dof_map,
    std::vector<ElementDofMap> &dof_map_reverse,
    std::vector<libMesh::dof_id_type> &handle_to_element_id,
    std::vector<ElementDofMap> &variable_order, 
    const libMesh::Parallel::Communicator &comm)
{
    using namespace std;

    //READ 1
    int dim, role;
    is >> dim >> role;

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
    handle_to_element_id.resize(n_elements);

    dof_map_reverse.resize(n_elements);

    CHECK_STREAM_READ_BEGIN("elements", is);

    for(long i = 0; i !=n_elements; ++i) {
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

        //READ 9
        is >> dof_map.at(i);
        is >> dof_map_reverse.at(i);

        mesh_ptr->add_elem(elem);

        libmesh_assert(elem);

    }

    CHECK_STREAM_READ_END("elements", is);

    //READ 10
    //  variable_number.resize(1);
    //  is >> variable_number.at(0);

    //READ 11
    variable_order.resize(1);
    is >> variable_order.at(0);

    //!!!! dummy parameters
    space = mesh_ptr;
}

static void read_spaces(
    moonolith::InputStream &is,
    FESpacesRAdapter &spaces,
    const libMesh::Parallel::Communicator &comm_master,
    const libMesh::Parallel::Communicator &comm_slave)
{

    bool has_master, has_slave;
    is >> has_master >> has_slave;

    spaces.spaces().resize(2);


    if(has_master) {
        read_space(
            is,
            spaces.spaces()[0],
            spaces.dof_map(0),
            spaces.dof_map_reverse(0),
            spaces.handle_to_element_id(0), 
            spaces.variable_order(0),
            comm_master);

        } else {
            spaces.spaces()[0] = nullptr;
        }

        if(has_slave) {
            read_space(
                is, 
                spaces.spaces()[1],
                spaces.dof_map(1),
                spaces.dof_map_reverse(1),
                spaces.handle_to_element_id(1),
                spaces.variable_order(1),
                comm_slave);

            } else {
                spaces.spaces()[1] = nullptr;
            }
        }
    }

#endif //UTOPIA_FE_SPACES_R_ADAPTER_HPP
