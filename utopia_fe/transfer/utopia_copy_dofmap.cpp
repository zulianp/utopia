
#include "utopia_copy_dofmap.hpp"

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

  Copy_Dof_Map::Copy_Dof_Map(){}

    void  Copy_Dof_Map::copy_global_dofs(MeshBase &space,
                          const std::shared_ptr<libMesh::DofMap>  &original_dof_map,
                          const unsigned int  &var_num,
                          std::vector<ElementDofMap> &dof_map,
                          std::vector<ElementDofMap> &variable_type, 
                          const int n_elements){
            

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
        
      
      void   Copy_Dof_Map::copy_var_order(DofMap &dofmap, std::vector<ElementDofMap> &variable_order)
      {
            variable_order.resize(1);
            FEType fe_type =  dofmap.variable(0).type();
            variable_order[0].global.push_back(fe_type.order);
      }
    }