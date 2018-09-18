#ifndef LibmeshContactForMoose_HPP
#define LibmeshContactForMoose_HPP

#include "utopia_assemble_contact.hpp"
#include "moonolith_communicator.hpp"

#include <iostream>

namespace utopia {

   inline bool MooseSurfaceAssemble(moonolith::Communicator &comm,
    const std::shared_ptr<libMesh::MeshBase> &mesh,
    const std::shared_ptr<libMesh::DofMap> &dof_map,
    const std::shared_ptr<const unsigned int> &_var_num,
    USparseMatrix &B,
    USparseMatrix &orthogonal_trafos,
    UVector &gap,
    USparseMatrix &normals,
    UVector &is_contact_node,
    const libMesh::Real search_radius,
    const std::vector< std::pair<int, int> > &tags,
    const bool use_biorth = true)
   {
       return assemble_contact(comm, mesh, dof_map, *_var_num, B, orthogonal_trafos, gap, normals, is_contact_node, search_radius, tags, use_biorth);
   }


   inline bool MooseSurfaceAssemble(
    moonolith::Communicator &comm,
    const std::shared_ptr<libMesh::MeshBase> &mesh,
    const std::shared_ptr<libMesh::DofMap> &dof_map,
    const std::shared_ptr<const unsigned int> &_var_num,
    USparseMatrix &B,
    USparseMatrix &orthogonal_trafos,
    UVector &gap,
    USparseMatrix &normals,
    UVector &is_contact_node,
    const libMesh::Real search_radius,
    const int tag_1, 
    const int tag_2,
    const bool use_biorth = true)
    {
        return assemble_contact(comm, mesh, dof_map, *_var_num, B, orthogonal_trafos, gap, normals, is_contact_node, search_radius, tag_1, tag_2, use_biorth);
    }
}



#endif //LIBMESH_CUTLIBPP_ADAPTERS_HPP



