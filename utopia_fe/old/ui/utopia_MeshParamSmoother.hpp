#ifndef UTOPIA_MESH_PARAM_SMOOTHER_HPP
#define UTOPIA_MESH_PARAM_SMOOTHER_HPP

#include "utopia_libmesh.hpp"
#include "utopia_ui.hpp"

namespace utopia {

    class MeshParamSmoother final : public Configurable {
    public:
        void read(Input &is) override;
        void apply(libMesh::MeshBase &mesh);

    private:
        int operator_power_ = 1;
        std::set<libMesh::boundary_id_type> requested_boundary_ids_;

        void apply_aux(const libMesh::MeshBase &surf_mesh,
                       const libMesh::DofMap &surf_dof_map,
                       const USparseMatrix &permutation,
                       const libMesh::DofMap &dof_map,
                       libMesh::MeshBase &mesh);
    };

}  // namespace utopia

#endif  // UTOPIA_MESH_PARAM_SMOOTHER_HPP
