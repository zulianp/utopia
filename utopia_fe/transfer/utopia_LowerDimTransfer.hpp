#ifndef UTOPIA_LOWER_DIM_TRANSFER_HPP
#define UTOPIA_LOWER_DIM_TRANSFER_HPP

#include "utopia_NewTransferAssembler.hpp"

namespace utopia {

    class LowerDimTransfer {
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
                // Q(std::make_shared<USparseMatrix>()),
                T(std::make_shared<USparseMatrix>()) {}

            void permute(const USparseMatrix &P, TransferData &out);

            std::shared_ptr<USparseMatrix> B, D, /*Q,*/ T;
        };

        bool assemble(
            const MeshBase &mesh,
            const DofMap   &dofs,
            const TransferOptions &opts = TransferOptions()
        );

        inline std::shared_ptr<PseudoL2TransferOperator> build_operator() const
        {
            return std::make_shared<PseudoL2TransferOperator>(data.T);
        }

        LowerDimTransfer() {}

    private:
        TransferData data;
    };
}

#endif //UTOPIA_LOWER_DIM_TRANSFER_HPP
