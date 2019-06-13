#ifndef UTOPIA_LIBMESH_DOF_MAP_ADAPTER_NEW_HPP
#define UTOPIA_LIBMESH_DOF_MAP_ADAPTER_NEW_HPP

#include "utopia_ElementDofMap.hpp"
#include "utopia_fe_base.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"
#include "moonolith_communicator.hpp"
#include "moonolith_function_space.hpp"

#include <vector>
#include <memory>

#include "libmesh/mesh.h"
#include "libmesh/elem.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/boundary_mesh.h"
#include "libmesh/dof_map.h"

namespace utopia {

    //TODO FEType
    class ElementDofMapAdapter {
    public:
        void init(
            libMesh::MeshBase &vol_mesh,
            const libMesh::DofMap &vol_dof_map,
            libMesh::MeshBase &surf_mesh,
            const libMesh::DofMap &surf_dof_map,
            const int var_num,
            const int surf_var_num);

    private:
        moonolith::DofMap dof_map_;
    };

}

#endif //UTOPIA_LIBMESH_DOF_MAP_ADAPTER_NEW_HPP
