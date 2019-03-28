#include "utopia_assemble_volume_transfer.hpp"
#include "utopia_L2LocalAssembler.hpp"
#include "utopia_InterpolationLocalAssembler.hpp"
#include "utopia_Local2Global.hpp"
#include "utopia_TransferAssembler.hpp"

#include <cmath>
#include <queue>
#include <algorithm>
#include <sstream>
#include <numeric>

using namespace libMesh;

namespace utopia {

    bool assemble_volume_transfer(
        moonolith::Communicator &,
        const std::shared_ptr<MeshBase> &from_mesh,
        const std::shared_ptr<MeshBase> &to_mesh,
        const std::shared_ptr<DofMap>   &from_dofs,
        const std::shared_ptr<DofMap>   &to_dofs,
        const unsigned int              &from_var_num,
        const unsigned int              &to_var_num,
        bool use_biorth,
        int n_var,
        USparseMatrix &B,
        const std::vector< std::pair<int, int> > &tags,
        const bool use_interpolation)
    {
        std::shared_ptr<LocalAssembler> assembler;

        if(use_interpolation) {
            std::cout << "[Status] using interpolation" << std::endl;
            assembler = std::make_shared<InterpolationLocalAssembler>(from_mesh->mesh_dimension());
        } else {
            std::cout << "[Status] using projection" << std::endl;
            assembler = std::make_shared<L2LocalAssembler>(from_mesh->mesh_dimension(), use_biorth);
        }

        auto local2global = std::make_shared<Local2Global>(use_interpolation);

        TransferOptions opts;
        opts.from_var_num = from_var_num;
        opts.to_var_num   = to_var_num;
        opts.n_var        = n_var;
        opts.tags         = tags;

        //////////////////////////////////////////////////////

        TransferAssembler transfer_assembler(assembler, local2global);
        return transfer_assembler.assemble(from_mesh, from_dofs, to_mesh, to_dofs, B, opts);
    }

}
