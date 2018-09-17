#ifndef UTOPIA_ASSEMBLE_VOLUME_TRANSFER_HPP
#define UTOPIA_ASSEMBLE_VOLUME_TRANSFER_HPP

#include <vector>
#include <utility>
#include <memory>

#include "utopia.hpp"
#include "utopia_fe_base.hpp"
#include "libmesh/libmesh_common.h"

//forward decl
namespace moonolith {
	class Communicator;
}

namespace libMesh {
	class MeshBase;
	class DofMap;
}

namespace utopia {

  bool assemble_volume_transfer(
      moonolith::Communicator &comm,
      const std::shared_ptr<libMesh::MeshBase> &master,
      const std::shared_ptr<libMesh::MeshBase> &slave,
      const std::shared_ptr<libMesh::DofMap> &dof_master,
      const std::shared_ptr<libMesh::DofMap> &dof_slave,
      const unsigned int &from_var_num,
      const unsigned int &to_var_num,
      bool use_biorth,
      int n_var,
      USMatrix &B,
      const std::vector< std::pair<int, int> > &tags = std::vector< std::pair<int, int> >(),
      const bool use_interpolation = false);
}

#endif //UTOPIA_ASSEMBLE_VOLUME_TRANSFER_HPP
