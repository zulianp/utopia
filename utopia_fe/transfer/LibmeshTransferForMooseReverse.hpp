
#ifndef LibmeshTransferForMoose_REVERSE_HPP
#define LibmeshTransferForMoose_REVERSE_HPP

#include "utopia_assemble_volume_transfer_r.hpp"
#include "moonolith_communicator.hpp"

namespace utopia {
        
    
    inline bool AssembleMOOSEReverse(moonolith::Communicator &comm,
                                    const std::shared_ptr<libMesh::MeshBase> &master,
                                    const std::shared_ptr<libMesh::MeshBase> &slave,
                                    const std::shared_ptr<libMesh::DofMap> &dof_master,
                                    const std::shared_ptr<libMesh::DofMap> &dof_slave,
                                    const std::shared_ptr<libMesh::DofMap> &dof_reverse_master,
                                    const std::shared_ptr<libMesh::DofMap> &dof_reverse_slave,
                                    const unsigned int & _from_var_num,
                                    const unsigned int & _to_var_num,
                                    const unsigned int & _from_var_num_r,
                                    const unsigned int & _to_var_num_r,
                                    std::unordered_set<utopia::SizeType> &particip_slave_dofs,
                                    bool  use_biorth_,
                                    int n_var,
                                    int n_var_r,
                                    DSMatrixd &B,
                                    DSMatrixd &B_reverse)
    {

      return assemble_volume_transfer_r(comm, 
                                        master, 
                                        slave, 
                                        dof_master, 
                                        dof_slave, 
                                        dof_reverse_master, 
                                        dof_reverse_slave, 
                                        _from_var_num, 
                                        _to_var_num, 
                                        _from_var_num_r, 
                                        _to_var_num_r,  
                                        particip_slave_dofs,
                                        use_biorth_, 
                                        n_var, 
                                        n_var_r, 
                                        B, 
                                        B_reverse);
    }

}


#endif //LIBMESH_TRANSFER_FOR_MOOSE_HPP
