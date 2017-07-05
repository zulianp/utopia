
#ifndef LibmeshTransferForMoose_HPP
#define LibmeshTransferForMoose_HPP

#include "utopia_assemble_volume_transfer.hpp"


namespace utopia {
        
    
    inline bool AssembleMOOSE(express::Communicator &comm,
                              const std::shared_ptr<libMesh::MeshBase> &master,
                              const std::shared_ptr<libMesh::MeshBase> &slave,
                              const std::shared_ptr<libMesh::DofMap> &dof_master,
                              const std::shared_ptr<libMesh::DofMap> &dof_slave,
                              const unsigned int & _from_var_num,
                              const unsigned int & _to_var_num,
                              bool  use_biorth_, int n_var, DSMatrixd &B)
    {
      return assemble_volume_transfer(comm, master, slave, dof_master, dof_slave, _from_var_num, _to_var_num, use_biorth_, n_var, B);
    }
}


#endif //LIBMESH_TRANSFER_FOR_MOOSE_HPP


