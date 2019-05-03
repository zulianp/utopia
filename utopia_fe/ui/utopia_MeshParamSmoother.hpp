#ifndef UTOPIA_MESH_PARAM_SMOOTHER_HPP
#define UTOPIA_MESH_PARAM_SMOOTHER_HPP

#include "utopia_ui.hpp"
#include "utopia_libmesh.hpp"

namespace utopia {

    class MeshParamSmoother final : public Configurable   {
    public:

        void read(Input &is) override;
        void apply(libMesh::MeshBase &mesh);

    private:
        int operator_power_ = 1;
    };

}

#endif //UTOPIA_MESH_PARAM_SMOOTHER_HPP
