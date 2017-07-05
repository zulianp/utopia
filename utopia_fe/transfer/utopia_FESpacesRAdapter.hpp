#ifndef UTOPIA_FE_SPACES_R_ADAPTER_HPP
#define UTOPIA_FE_SPACES_R_ADAPTER_HPP 

#include "utopia_ElementDofMap.hpp"
#include "Array.hpp"
#include "express_Communicator.hpp"

#include "libmesh/serial_mesh.h"
#include "libmesh/dof_map.h"

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

		inline FESpacesRAdapter(const express::Communicator &comm) : comm(comm){}
        	       
		FESpacesRAdapter(const std::shared_ptr<libMesh::MeshBase> &master,
                         const std::shared_ptr<libMesh::MeshBase> &slave,
		                 const std::shared_ptr<libMesh::DofMap>  &dof_map_master,
		                 const std::shared_ptr<libMesh::DofMap>  &dof_map_slave,
		                 const std::shared_ptr<libMesh::DofMap> &dof_map_reverse_master,
		                 const std::shared_ptr<libMesh::DofMap> &dof_map_reverse_slave,
		                 const unsigned int &_from_var_num,
		                 const unsigned int &_to_var_num,
		                 const unsigned int &_from_var_num_r,
		                 const unsigned int &_to_var_num_r);
        
        
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
        
        inline void set_must_destroy_attached(const int index, const bool value)
        {
            assert(index < 2);
            assert(index >= 0);
            must_destroy_attached[index] = value;
        }
        
        
        inline  std::vector<ElementDofMap> &variable_number(const int i)
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
        
        
    private:
        express::Communicator comm;
        std::vector<std::shared_ptr< libMesh::MeshBase>> spaces_;
        std::vector<ElementDofMap> dof_maps_[2];
        std::vector<ElementDofMap> dof_maps_reverse_[2];
        std::vector<ElementDofMap> var_number_[2];
        std::vector<ElementDofMap> var_order_[2];
        std::vector<ElementDofMap> var_type_[2];
        bool must_destroy_attached[2];
        
        
        
        
        static void copy_global_dofs_r(libMesh::MeshBase &space,
                                       const std::shared_ptr<libMesh::DofMap>  &original_dof_map,
                                       const std::shared_ptr<libMesh::DofMap>  &original_dof_map_reverse,
                                       const unsigned int  &var_num,
                                       const unsigned int  &var_num_r,
                                       std::vector<ElementDofMap> &dof_map,
                                       std::vector<ElementDofMap> &dof_map_reverse,
                                       std::vector<ElementDofMap> &variable_type, const int n_elements);
        
        static void copy_var_order_r(libMesh::DofMap &dofmap, std::vector<ElementDofMap> &variable_order);
        
    };
}

#endif //UTOPIA_FE_SPACES_R_ADAPTER_HPP
    