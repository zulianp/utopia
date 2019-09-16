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
            const libMesh::MeshBase &mesh,
            const libMesh::DofMap &dof_map,
            const int var_num);

        void init_surf_to_vol(
            libMesh::MeshBase &vol_mesh,
            const libMesh::DofMap &vol_dof_map,
            libMesh::MeshBase &surf_mesh,
            const libMesh::DofMap &surf_dof_map,
            const int var_num,
            const int surf_var_num);

        void make_element_node_map(ElementDofMapAdapter &out) const;

        static void make_permutation(const ElementDofMapAdapter &from, const ElementDofMapAdapter &to, USparseMatrix &mat);
        
        //original dof structure has to have the vector structure already
        //from vector to vector. from map has a scalar map
        static void make_vector_permutation(
            const int dim,
            const ElementDofMapAdapter &from,
            const ElementDofMapAdapter &to,
             USparseMatrix &mat);

        //original dof structure has to have the vector structure already
        //from scalar to repeated vector
        static void make_tensorize_permutation(
            const int dim,
            const ElementDofMapAdapter &from,
            const ElementDofMapAdapter &to,
             USparseMatrix &mat);

        inline const moonolith::Storage<moonolith::Dofs> &dofs() const {
            return dof_map_.dofs();
        }

        inline moonolith::Storage<moonolith::Dofs> &dofs() {
            return dof_map_.dofs();
        }

        inline SizeType n_local_dofs() const { return n_local_dofs_;}

    private:
        moonolith::Communicator comm_;
        SizeType n_local_dofs_;
        SizeType max_nnz_;
        moonolith::DofMap dof_map_;
    };

}

#endif //UTOPIA_LIBMESH_DOF_MAP_ADAPTER_NEW_HPP
