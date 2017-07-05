#include "utopia_FESpacesAdapter.hpp"

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

   FESpacesAdapter::FESpacesAdapter(
       const std::shared_ptr<MeshBase> &master,
       const std::shared_ptr<MeshBase> &slave,
       const std::shared_ptr<libMesh::DofMap>  &dof_map_master,
       const std::shared_ptr<libMesh::DofMap>  &dof_map_slave,
       const unsigned int &_from_var_num,
       const unsigned int &_to_var_num)
   {
               
               spaces_.reserve(2);
               spaces_.push_back(master);
               spaces_.push_back(slave);

               const int n_elements_master = master->n_elem();
               const int n_elements_slave  = slave->n_elem();
               
               const int n_elements = n_elements_master + n_elements_slave;

               
               std::cout<<"MASTER DOF"<<std::endl;
               copy_global_dofs(*master, dof_map_master, _from_var_num, dof_maps_[0],  var_type_[0], n_elements);
               
               std::cout<<"SLAVE DOF"<<std::endl;
               copy_global_dofs(*slave, dof_map_slave,  _to_var_num,  dof_maps_[1], var_type_[1], n_elements);

               copy_var_order(*dof_map_master, var_order_[0]);
               copy_var_order(*dof_map_slave,  var_order_[1]);
               
      }

    void copy_global_dofs(MeshBase &space,
                          const std::shared_ptr<libMesh::DofMap>  &original_dof_map,
                          const unsigned int  &var_num,
                          std::vector<ElementDofMap> &dof_map,
                          std::vector<ElementDofMap> &variable_type, 
                          const int n_elements)
    {
            

            std::vector<dof_id_type> temp;
 
            dof_map.resize(n_elements);
            
            MeshBase::const_element_iterator e_it        =space.active_local_elements_begin();
            const MeshBase::const_element_iterator e_end =space.active_local_elements_end();
            
            int i=0;
            
            variable_type.resize(1);
            
            bool first=true;
            
            
            for (; e_it != e_end; ++e_it){
                
                Elem *elem = *e_it;
                
                original_dof_map->dof_indices(elem, temp, var_num);
                
                dof_map[elem->id()].global.insert(dof_map[elem->id()].global.end(), temp.begin(), temp.end());
                
                
                if (first)
                {
                    variable_type[0].global.push_back(elem->type());
                    first=false;
                    
                }
                i++;
            }
            
        }
        
      
      void copy_var_order(DofMap &dofmap, std::vector<ElementDofMap> &variable_order)
      {
            variable_order.resize(1);
            FEType fe_type =  dofmap.variable(0).type();
            variable_order[0].global.push_back(fe_type.order);
      }
}