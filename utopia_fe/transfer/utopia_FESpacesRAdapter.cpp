#include "utopia_FESpacesRAdapter.hpp"

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

   FESpacesRAdapter::FESpacesRAdapter(
       const std::shared_ptr<MeshBase> &master,
       const std::shared_ptr<MeshBase> &slave,
       const std::shared_ptr<libMesh::DofMap>  &dof_map_master,
       const std::shared_ptr<libMesh::DofMap>  &dof_map_slave,
       const std::shared_ptr<DofMap> &dof_map_reverse_master,
       const std::shared_ptr<DofMap> &dof_map_reverse_slave,
       const unsigned int &_from_var_num,
       const unsigned int &_to_var_num,
       const unsigned int &_from_var_num_r,
       const unsigned int &_to_var_num_r)
   {
               
               spaces_.reserve(2);
               spaces_.push_back(master);
               spaces_.push_back(slave);

               const int n_elements_master = master->n_elem();
               const int n_elements_slave  = slave->n_elem();
               
               const int n_elements = n_elements_master + n_elements_slave;
               
               std::cout<<"MASTER DOF"<<std::endl;
               copy_global_dofs_r( *master ,dof_map_master, dof_map_reverse_master, _from_var_num, _to_var_num_r, 
                                   dof_maps_[0],  dof_maps_reverse_[0], var_type_[0], n_elements);
               
               std::cout<<"SLAVE DOF"<<std::endl;
               copy_global_dofs_r( *slave,dof_map_slave,   dof_map_reverse_slave,  _to_var_num,   _from_var_num_r,   
                                   dof_maps_[1], dof_maps_reverse_[1], var_type_[1], n_elements);
          
               copy_var_order_r(*dof_map_master, var_order_[0]);
               copy_var_order_r(*dof_map_slave,  var_order_[1]);
               
      }

       void  FESpacesRAdapter::copy_global_dofs_r(MeshBase &space,
                             const std::shared_ptr<DofMap>  &original_dof_map,
                             const std::shared_ptr<DofMap>  &original_dof_map_reverse,
                             const unsigned int  &var_num,
                             const unsigned int  &var_num_r,
                             std::vector<ElementDofMap> &dof_map,
                             std::vector<ElementDofMap> &dof_map_reverse,
                             std::vector<ElementDofMap> &variable_type, const int n_elements)
       {
               
               //            auto &mesh = space.get_mesh();
               //            auto &original_dof_map = space.get_dof_map();
               std::vector<dof_id_type> temp;
               std::vector<dof_id_type> temp_reverse;
               dof_map.resize(n_elements);
               dof_map_reverse.resize(n_elements);
               
               
               //        std::cout<<"______________________________COPY_DOF_BEGIN____________________________"<<std::endl;
               
               MeshBase::const_element_iterator e_it        =space.active_local_elements_begin();
               const MeshBase::const_element_iterator e_end =space.active_local_elements_end();
               
               int i=0;
               
               variable_type.resize(1);
               
               bool first=true;
               
               
               for (; e_it != e_end; ++e_it){
                   
                   Elem *elem = *e_it;
                   
                   original_dof_map->dof_indices(elem, temp, var_num);
                   
                   original_dof_map_reverse->dof_indices(elem, temp_reverse, var_num_r);
                   
                   dof_map[elem->id()].global.insert(dof_map[elem->id()].global.end(), temp.begin(), temp.end());
                   
                   dof_map_reverse[elem->id()].global.insert(dof_map_reverse[elem->id()].global.end(), temp_reverse.begin(), temp_reverse.end());
                   
                   if (first)
                   {
                       variable_type[0].global.push_back(elem->type());
                       first=false;
                       
                   }
                   i++;
               }
               
           }
        
      
      void   FESpacesRAdapter::copy_var_order_r(DofMap &dofmap, std::vector<ElementDofMap> &variable_order)
      {
            variable_order.resize(1);
            FEType fe_type =  dofmap.variable(0).type();
            variable_order[0].global.push_back(fe_type.order);
      }
 }