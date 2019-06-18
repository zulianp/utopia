#ifndef UTOPIA_NEW_TRANSFER_ASSEMBLER_HPP
#define UTOPIA_NEW_TRANSFER_ASSEMBLER_HPP

#include "utopia_TransferAssembler.hpp"
#include <memory>

namespace utopia {

    class TransferData {
    public:
        TransferData() : 
            B(std::make_shared<USparseMatrix>()),
            D(std::make_shared<USparseMatrix>()),
            Q(std::make_shared<USparseMatrix>()),
            T(std::make_shared<USparseMatrix>())
        {}

        std::shared_ptr<USparseMatrix> B, D, Q, T;
    };

    class NewTransferAssembler {
    public:
        using FunctionSpace = utopia::LibMeshFunctionSpace;
        using SparseMatrix  = utopia::USparseMatrix;
        using MeshBase      = libMesh::MeshBase;
        using DofMap        = libMesh::DofMap;

        bool assemble(
            const std::shared_ptr<MeshBase> &from_mesh,
            const std::shared_ptr<DofMap>   &from_dofs,
            const std::shared_ptr<MeshBase> &to_mesh,
            const std::shared_ptr<DofMap>   &to_dofs,
            const TransferOptions &opts = TransferOptions()
        );

        inline std::shared_ptr<PseudoL2TransferOperator> build_operator() const
        {
            return std::make_shared<PseudoL2TransferOperator>(data.T);
        }

        TransferData data;
    };
}

#endif //UTOPIA_NEW_TRANSFER_ASSEMBLER_HPP
