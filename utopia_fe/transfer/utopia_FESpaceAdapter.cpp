#include "utopia_FESpaceAdapter.hpp"

#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"

#include <numeric>
#include <cmath>

using libMesh::MeshBase;
using libMesh::Elem;
using libMesh::dof_id_type;
using libMesh::Node;
using libMesh::DofMap;
using libMesh::FEType;

namespace utopia {
    FESpaceAdapter::FESpaceAdapter(
        const std::shared_ptr<MeshBase> &master_slave,
        const std::shared_ptr<libMesh::DofMap>  &original_dofmap,
        const unsigned int var_num,
        const std::vector< std::pair<int, int>> &tags)
    {

        mesh_ = master_slave;

        moonolith::Communicator comm = master_slave->comm().get();

        copy_global_dofs(*master_slave,
                         original_dofmap,
                         var_num,
                         dof_maps_,
                         var_type_,
                         var_number_,
                         subdomain_id_,
                         side_set_id_,
                         side_set_id_tag_,
                         face_set_id_global_,
                         ownershipRangesFaceID_,
                         handle_to_element_id_,
                         tags);

        copy_var_order(*original_dofmap, var_order_);
    }

    bool FESpaceAdapter::is_tagged_contact_boundary(const MeshBase &mesh,
                                             const Elem *elem,
                                             const int side_elem,
                                             const std::vector< std::pair<int, int> > &tags)
    {
        for(auto t : tags) {
            if (mesh.get_boundary_info().has_boundary_id(elem, side_elem, t.first) ||
                mesh.get_boundary_info().has_boundary_id(elem, side_elem, t.second)) {
                 return true;
             }
         }

         return false;
     }

     void FESpaceAdapter::copy_global_dofs(MeshBase &mesh,
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
                                         const std::vector< std::pair<int, int> >  &tags)
     {
         std::vector<dof_id_type> temp;

         moonolith::Communicator comm = mesh.comm().get();

         ownershipRangesFaceID.resize(comm.size() + 1, 0);

         libMesh::dof_id_type n_elements = mesh.n_active_local_elem();
         libMesh::dof_id_type local_element_id = 0;

         std::vector<ElementDofMap> face_set_id;
         dof_map.resize(n_elements);
         subdomain_id.resize(n_elements);
         side_set_id.resize(n_elements);
         side_set_id_tag.resize(n_elements); //non serve
         face_set_id.resize(n_elements);
         face_set_id_global.resize(n_elements);
         handle_to_element_id.resize(n_elements);

         bool first=true;

         std::vector<const Node *> elem_nodes;
         int jj_side_id_one = 0;
         int jj_side_id_one_tag = 0;
         int jj_side_id_one_check = 0;
         int offset=0;
         int f_id=0;
//	     int n_f=0;
         MeshBase::const_element_iterator e_it_s = mesh.active_local_elements_begin();
         const MeshBase::const_element_iterator e_end_s = mesh.active_local_elements_end();

         for (; e_it_s != e_end_s; ++e_it_s, ++local_element_id) {
             Elem * elem = *e_it_s;

             handle_to_element_id[local_element_id] = elem->id();

             bool  check_side_id_one=true;
             bool  check_side_id_one_tag=true;
             bool  check_side_id_one_check=true;
//	         bool  check_face_id=true;

             for (int side_elem=0; side_elem<elem->n_sides(); side_elem++) {
                 if (check_side_id_one) {
                     check_side_id_one=false;
                     jj_side_id_one++;
                 }
             }

             if (elem->on_boundary()){
                 for (int side_elem=0; side_elem<elem->n_sides(); side_elem++) {
                     {
                         if (is_tagged_contact_boundary(mesh, elem, side_elem, tags) && check_side_id_one_tag){
                             side_set_id[local_element_id].global.push_back(mesh.get_boundary_info().boundary_id(elem, side_elem));
                             check_side_id_one_tag = false;
                             jj_side_id_one_tag++;
                         }
                     }
                 }
             }

             for (int side_elem=0; side_elem<elem->n_sides(); side_elem++) {
                 if (check_side_id_one_check){
                     check_side_id_one_check=false;
                     jj_side_id_one_check++;
                 }
             }

             if (elem->on_boundary()) {

                 for(unsigned int side_elem = 0; side_elem < elem->n_sides(); ++side_elem) {
                     {
                         if (is_tagged_contact_boundary(mesh, elem, side_elem, tags)) {

                             face_set_id[local_element_id].global.push_back(f_id++);

                             offset++;
                         } else {
                             face_set_id[local_element_id].global.push_back(-1);
                         }
                     }
                 }
             } else {
                 face_set_id[local_element_id].global.insert(face_set_id[local_element_id].global.end(), -1);
             }


             subdomain_id[local_element_id].global.insert(subdomain_id[local_element_id].global.end(),elem->subdomain_id());

             original_dofmap->dof_indices(elem, temp, var_num);

             dof_map[local_element_id].global.insert(dof_map[local_element_id].global.end(), temp.begin(), temp.end());

             if(first) {
                 //works only because libmesh does not support mixed elements
                 variable_type.push_back(elem->type());
                 variable_number.push_back(var_num);
                 first = false;
             }
         }

         ownershipRangesFaceID[comm.rank()+1]+= static_cast<unsigned int>(offset);

         comm.all_reduce(&ownershipRangesFaceID[0],  ownershipRangesFaceID.size(),  moonolith::MPISum());

         std::partial_sum(ownershipRangesFaceID.begin(), ownershipRangesFaceID.end(),
                          ownershipRangesFaceID.begin());

         MeshBase::const_element_iterator e_it_new = mesh.active_local_elements_begin();
         const MeshBase::const_element_iterator e_end_new = mesh.active_local_elements_end();

         local_element_id = 0;
         for (; e_it_new != e_end_new; ++e_it_new, ++local_element_id)
         {
             Elem * elem_new = *e_it_new;

             if (elem_new->on_boundary()){
                 for (int jj=0; jj<face_set_id[local_element_id].global.size(); jj++){
                     int i = face_set_id[local_element_id].global.at(jj);

                     if(i != -1) {
                         int global_id = i + ownershipRangesFaceID[comm.rank()];
                         face_set_id_global[local_element_id].global.insert(face_set_id_global[local_element_id].global.end(),global_id);
                     } else {
                         face_set_id_global[local_element_id].global.insert(face_set_id_global[local_element_id].global.end(),-1);
                     }
                 }
             }
         }
     }

     void FESpaceAdapter::copy_var_order(DofMap &dofmap, std::vector<long> &variable_order)
     {
         FEType fe_type =  dofmap.variable(0).type();
         variable_order.push_back(fe_type.order);
     }
}
