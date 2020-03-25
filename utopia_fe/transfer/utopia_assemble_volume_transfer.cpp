#include "utopia_assemble_volume_transfer.hpp"
#include "utopia_L2LocalAssembler.hpp"
#include "utopia_InterpolationLocalAssembler.hpp"
#include "utopia_Local2Global.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_MeshTransferOperator.hpp"
#include "utopia_ui.hpp"
#include "utopia_NewTransferAssembler.hpp"

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
        // from_mesh->comm().barrier();
        // if(from_mesh->comm().rank() == 0) { moonolith::logger() << "assemble_volume_transfer" << std::endl; }

#ifdef WITH_NEW_TRANSFER
#warning "enabling new transfer features with WITH_NEW_TRANSFER"
        if(use_biorth && !use_interpolation) {

            from_mesh->comm().barrier();
            if(from_mesh->comm().rank() == 0) { moonolith::logger() << "[Status] using new transfer" << std::endl; }

            TransferOptions opts;
            opts.from_var_num = from_var_num;
            opts.to_var_num   = to_var_num;
            opts.n_var        = n_var;
            opts.tags         = tags;

            NewTransferAssembler new_assembler;
            // new_assembler.use_convert_transfer(params_->use_convert_transfer);
            // new_assembler.remove_incomplete_intersections(params_->remove_incomplete_intersections);
            // new_assembler.handle_adaptive_refinement(params_->handle_adaptivity);
            
            // if(params_->surface_transfer) {
            //     if(!new_assembler.surface_assemble(from_mesh, from_dofs, to_mesh, to_dofs, opts)) {
            //         return false;
            //     }
            // } else {
                if(!new_assembler.assemble(from_mesh, from_dofs, to_mesh, to_dofs, opts)) {
                    return false;
                }
            // }

            auto op = new_assembler.build_operator();
            B = std::move(*op->matrix());
            return true;

        } 
#endif //WITH_NEW_TRANSFER

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
