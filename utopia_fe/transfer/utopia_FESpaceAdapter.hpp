#ifndef UTOPIA_FE_SPACE_ADAPTER_HPP
#define UTOPIA_FE_SPACE_ADAPTER_HPP 

#include "utopia_ElementDofMap.hpp"
#include "moonolith_communicator.hpp"

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

	private:
		moonolith::Communicator comm;
		std::shared_ptr< libMesh::MeshBase> mesh_;
		std::vector<ElementDofMap> dof_maps_;
		
		std::vector<long> var_number_;
		std::vector<long> var_type_;
		std::vector<long> var_order_;


		std::vector<ElementDofMap> subdomain_id_;

		std::vector<ElementDofMap> side_set_id_;
		std::vector<ElementDofMap> face_set_id_global_;
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
			const int n_elements, 
			std::vector<long> &variable_number,
			std::vector<ElementDofMap> &subdomain_id, 
			std::vector<ElementDofMap> &side_set_id,
			std::vector<ElementDofMap> &side_set_id_tag, 
			std::vector<ElementDofMap> &face_set_id_global,
			std::vector<moonolith::Integer> &ownershipRangesFaceID, 
			const std::vector< std::pair<int, int> > &tags);


		static void copy_var_order(libMesh::DofMap &dofmap, std::vector<long> &variable_order);
	};
}


#endif //UTOPIA_FE_SPACE_ADAPTER_HPP
