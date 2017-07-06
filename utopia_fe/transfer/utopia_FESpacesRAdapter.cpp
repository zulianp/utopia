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

               Copy_Dof_Map::copy_global_dofs( *slave, dof_map_slave,  _to_var_num,  
                                   dof_maps_[1], var_type_[1], n_elements);


               Copy_Dof_Map::copy_global_dofs( *master, dof_map_master, _from_var_num, 
                                   dof_maps_[0], var_type_[0], n_elements);
               

               Copy_Dof_Map::copy_global_dofs( *master , dof_map_reverse_master,  _to_var_num_r, 
                                 dof_maps_reverse_[0], var_type_[0], n_elements);
               

               Copy_Dof_Map::copy_global_dofs( *slave, dof_map_reverse_slave, _from_var_num_r,   
                                  dof_maps_reverse_[1], var_type_[1], n_elements);


          
               Copy_Dof_Map::copy_var_order(*dof_map_master, var_order_[0]);
               
               Copy_Dof_Map::copy_var_order(*dof_map_slave,  var_order_[1]);
               
      }

 }