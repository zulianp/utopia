#ifndef UTOPIA_NEW_TRANSFER_ASSEMBLER_HPP
#define UTOPIA_NEW_TRANSFER_ASSEMBLER_HPP

#include "utopia_TransferAssembler.hpp"
#include <memory>

namespace utopia {

    class NewTransferAssembler {
    public:
        using FunctionSpace = utopia::LibMeshFunctionSpace;
        using SparseMatrix  = utopia::USparseMatrix;
        using MeshBase      = libMesh::MeshBase;
        using DofMap        = libMesh::DofMap;

        class TransferData {
        public:
            TransferData() : 
                B(std::make_shared<USparseMatrix>()),
                D(std::make_shared<USparseMatrix>()),
                Q(std::make_shared<USparseMatrix>()),
                T(std::make_shared<USparseMatrix>())
            {}

            void permute(const USparseMatrix &P, TransferData &out);

            std::shared_ptr<USparseMatrix> B, D, Q, T;
        };

        bool assemble(
            const std::shared_ptr<MeshBase> &from_mesh,
            const std::shared_ptr<DofMap>   &from_dofs,
            const std::shared_ptr<MeshBase> &to_mesh,
            const std::shared_ptr<DofMap>   &to_dofs,
            const TransferOptions &opts = TransferOptions()
        );

        bool assemble(
            const MeshBase &from_mesh,
            const DofMap   &from_dofs,
            const MeshBase &to_mesh,
            const DofMap   &to_dofs,
            const TransferOptions &opts  = TransferOptions()
        );

        bool surface_assemble(
            const std::shared_ptr<MeshBase> &from_mesh,
            const std::shared_ptr<DofMap>   &from_dofs,
            const std::shared_ptr<MeshBase> &to_mesh,
            const std::shared_ptr<DofMap>   &to_dofs,
            const TransferOptions &opts = TransferOptions()
        );

        bool surface_assemble(
            const MeshBase &mesh,
            const DofMap   &dofs,
            const TransferOptions &opts = TransferOptions()
        );

        inline std::shared_ptr<PseudoL2TransferOperator> build_operator() const
        {
            return std::make_shared<PseudoL2TransferOperator>(data.T);
        }

        NewTransferAssembler() : use_convert_transfer_(true), remove_incomplete_intersections_(false) {}

        inline void use_convert_transfer(const bool val)
        {
            use_convert_transfer_ = val;
        }

        inline void remove_incomplete_intersections(const bool val)
        {
            remove_incomplete_intersections_ = val;
        }

    private:
        TransferData data;
        bool use_convert_transfer_;
        bool remove_incomplete_intersections_;
    };
}

#endif //UTOPIA_NEW_TRANSFER_ASSEMBLER_HPP
