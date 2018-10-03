#ifndef UTOPIA_ASSEMBLE_VOLUME_TRANSFER_R_HPP 
#define UTOPIA_ASSEMBLE_VOLUME_TRANSFER_R_HPP 

#include <vector>
#include <utility>
#include <memory>

#include <Eigen/Core>

#include "utopia.hpp"
#include "utopia_fe_base.hpp"
#include "libmesh/libmesh_common.h"

#include "moonolith_communicator.hpp"


namespace libMesh {
      class MeshBase;
      class DofMap;
}

namespace utopia {

      bool assemble_volume_transfer_r(
            moonolith::Communicator &comm,
            const std::shared_ptr<libMesh::MeshBase> &master,
            const std::shared_ptr<libMesh::MeshBase> &slave,
            const std::shared_ptr<libMesh::DofMap> &dof_master,
            const std::shared_ptr<libMesh::DofMap> &dof_slave,
            const std::shared_ptr<libMesh::DofMap> &dof_reverse_master,
            const std::shared_ptr<libMesh::DofMap> &dof_reverse_slave,
            const unsigned int &from_var_num,
            const unsigned int &to_var_num,
            const unsigned int &from_var_num_r,
            const unsigned int &to_var_num_r,
            bool use_biorth,
            int n_var,
            int n_var_r,
            USparseMatrix &B,
            USparseMatrix &B_reverse,
            const std::vector< std::pair<int, int> > &tags=std::vector< std::pair<int, int> >());

}

#endif //UTOPIA_ASSEMBLE_VOLUME_TRANSFER_R_HPP
