#ifndef UTOPIA_ASSEMBLE_VOLUME_TRANSFER_R_HPP 
#define UTOPIA_ASSEMBLE_VOLUME_TRANSFER_R_HPP 

#include <vector>
#include <utility>
#include <memory>

#include "utopia.hpp"
#include "libmesh/libmesh_common.h"

//forward decl
namespace express {
	class Communicator;
}

namespace libMesh {
	class MeshBase;
	class DofMap;
}

namespace utopia {


     bool assemble_volume_transfer_r(express::Communicator &comm,
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
                                    bool  use_biorth_,
                                    int n_var,
                                    int n_var_r,
                                    DSMatrixd &B,
                                    DSMatrixd &B_reverse);
	

}

#endif //UTOPIA_ASSEMBLE_VOLUME_TRANSFER_R_HPP 