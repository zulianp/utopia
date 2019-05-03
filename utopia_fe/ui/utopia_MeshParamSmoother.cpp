#include "utopia_MeshParamSmoother.hpp"

namespace utopia {

    void MeshParamSmoother::read(Input &is)
    {
        is.get("operator-power", operator_power_);

    }

    void MeshParamSmoother::apply(libMesh::MeshBase &mesh)
    {
        const auto spatial_dim = mesh.spatial_dimension();
        const auto mesh_dim = mesh.mesh_dimension();

        const bool is_surf = spatial_dim > mesh_dim;

        LibMeshFunctionSpace V(mesh, libMesh::LAGRANGE, libMesh::SECOND);
        V.initialize();

        auto &dof_map = V.dof_map();

        std::vector<bool> node_is_bdry;
        std::vector<UVector> rhs(spatial_dim);

        for(auto &r : rhs) {
          r = local_zeros(dof_map.n_local_dofs());
        }

        if(is_surf) {
            assert(false && "IMPLEMENT ME");
        } else {
            
            for(auto e_it = elements_begin(mesh); e_it != elements_end(mesh); ++e_it) {
                auto &e = **e_it;

                auto n_sides = e.n_sides();

                for(uint i = 0; i < n_sides; ++i) {
                    if(e.neighbor_ptr(i) == nullptr) {

                        auto side_ptr = e.build_side_ptr(i);

                    }
                }
            } 
        }

    }

}
