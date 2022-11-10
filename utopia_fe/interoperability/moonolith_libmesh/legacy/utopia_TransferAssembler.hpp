#ifndef UTOPIA_TRANSFER_ASSEMBLER_HPP
#define UTOPIA_TRANSFER_ASSEMBLER_HPP

#include "utopia.hpp"
#include "utopia_Local2Global.hpp"
#include "utopia_LocalAssembler.hpp"
#include "utopia_Path.hpp"
// #include "utopia_libmesh_old.hpp"
#include "utopia_ui.hpp"

#include "utopia_libmesh_TransferOperator.hpp"

#include <memory>

namespace utopia {

    class TransferAssembler final {
    public:
        using SparseMatrix = utopia::USparseMatrix;
        using MeshBase = libMesh::MeshBase;
        using DofMap = libMesh::DofMap;

        class Algorithm {
        public:
            virtual ~Algorithm() {}
            // virtual bool assemble(SparseMatrix &B) = 0;
            virtual bool assemble(std::vector<std::shared_ptr<SparseMatrix>> &B) = 0;
        };

        TransferAssembler(const std::shared_ptr<LocalAssembler> &assembler,
                          const std::shared_ptr<Local2Global> &local2global);

        bool assemble(const std::shared_ptr<MeshBase> &from_mesh,
                      const std::shared_ptr<DofMap> &from_dofs,
                      const std::shared_ptr<MeshBase> &to_mesh,
                      const std::shared_ptr<DofMap> &to_dofs,
                      SparseMatrix &B,
                      const TransferOptions &opts = TransferOptions());

        bool assemble(const std::shared_ptr<MeshBase> &from_mesh,
                      const std::shared_ptr<DofMap> &from_dofs,
                      const std::shared_ptr<MeshBase> &to_mesh,
                      const std::shared_ptr<DofMap> &to_dofs,
                      std::vector<std::shared_ptr<SparseMatrix>> &B,
                      const TransferOptions &opts = TransferOptions());

        void set_assembler(const std::shared_ptr<LocalAssembler> &assembler) { assembler_ = assembler; }

        void set_local_2_global(const std::shared_ptr<Local2Global> &local2global) { local2global_ = local2global; }

        ~TransferAssembler();

    private:
        std::shared_ptr<LocalAssembler> assembler_;
        std::shared_ptr<Local2Global> local2global_;
        std::shared_ptr<Algorithm> algorithm_;
    };

}  // namespace utopia

#endif  // UTOPIA_TRANSFER_ASSEMBLER_HPP
