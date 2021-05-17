#ifndef UTOPIA_GRID_2_MESH_TRANSFER_APP
#define UTOPIA_GRID_2_MESH_TRANSFER_APP

#include "utopia_FEApp.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_libmesh_old.hpp"

#include <memory>
#include <string>

namespace utopia {
    class LocalAssembler;
    class Local2Global;

    class Grid2MeshTransferApp final : public FEApp {
    public:
        class InputSpace;
        class InputGrid;

        ~Grid2MeshTransferApp();
        Grid2MeshTransferApp();

        void run(Input &in) override;

        static std::string command() { return "-g2m_transfer"; }
    };
}  // namespace utopia

#endif  // UTOPIA_GRID_2_MESH_TRANSFER_APP
