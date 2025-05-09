#ifndef UTOPIA_NEW_TRANSFER_ASSEMBLER_HPP
#define UTOPIA_NEW_TRANSFER_ASSEMBLER_HPP

#include <libmesh/dof_map.h>

#include <libmesh/mesh_base.h>
#include "utopia_fe_base.hpp"

#include <memory>

#include "utopia_libmesh_TransferOperator.hpp"

namespace utopia {

    class NewTransferAssembler {
    public:
        // using FunctionSpace = utopia::LibMeshFunctionSpace;
        using SparseMatrix = utopia::USparseMatrix;
        using MeshBase = libMesh::MeshBase;
        using DofMap = libMesh::DofMap;

        class TransferData {
        public:
            TransferData()
                : B(std::make_shared<USparseMatrix>()),
                  D(std::make_shared<USparseMatrix>()),
                  Q(std::make_shared<USparseMatrix>()),
                  T(std::make_shared<USparseMatrix>()),
                  constraint_matrix_from(std::make_shared<USparseMatrix>()),
                  constraint_matrix_to(std::make_shared<USparseMatrix>()),
                  post_constraint_matrix_from(std::make_shared<USparseMatrix>()),
                  post_constraint_matrix_to(std::make_shared<USparseMatrix>()) {}

            void permute(const USparseMatrix &P, TransferData &out);

            std::shared_ptr<USparseMatrix> B, D, Q, T;
            std::shared_ptr<USparseMatrix> constraint_matrix_from, constraint_matrix_to;
            std::shared_ptr<USparseMatrix> post_constraint_matrix_from, post_constraint_matrix_to;
        };

        bool assemble(const std::shared_ptr<MeshBase> &from_mesh,
                      const std::shared_ptr<DofMap> &from_dofs,
                      const std::shared_ptr<MeshBase> &to_mesh,
                      const std::shared_ptr<DofMap> &to_dofs,
                      const TransferOptions &opts = TransferOptions());

        bool assemble(const MeshBase &from_mesh,
                      const DofMap &from_dofs,
                      const MeshBase &to_mesh,
                      const DofMap &to_dofs,
                      const TransferOptions &opts = TransferOptions());

        bool surface_assemble(const std::shared_ptr<MeshBase> &from_mesh,
                              const std::shared_ptr<DofMap> &from_dofs,
                              const std::shared_ptr<MeshBase> &to_mesh,
                              const std::shared_ptr<DofMap> &to_dofs,
                              const TransferOptions &opts = TransferOptions());

        bool surface_assemble(const MeshBase &mesh,
                              const DofMap &dofs,
                              const TransferOptions &opts = TransferOptions());

        inline std::shared_ptr<PseudoL2TransferOperator> build_operator() const {
            return std::make_shared<PseudoL2TransferOperator>(data.T);
        }

        NewTransferAssembler()
            : use_convert_transfer_(true),
              remove_incomplete_intersections_(false),
              handle_adaptive_refinement_(false),
              use_dual_lagrange_multiplier_(true) {}

        inline void use_convert_transfer(const bool val) { use_convert_transfer_ = val; }

        inline void remove_incomplete_intersections(const bool val) { remove_incomplete_intersections_ = val; }

        inline void handle_adaptive_refinement(const bool val) { handle_adaptive_refinement_ = val; }

        inline void constraint_matrix_from(const std::shared_ptr<USparseMatrix> &mat) {
            data.constraint_matrix_from = mat;
        }

        inline void constraint_matrix_to(const std::shared_ptr<USparseMatrix> &mat) { data.constraint_matrix_to = mat; }

        inline void use_dual_lagrange_multiplier(const bool val) { use_dual_lagrange_multiplier_ = val; }

        std::vector<std::shared_ptr<USparseMatrix>> matrices() const {
            if (data.Q) {
                return {data.B, data.D, data.Q};
            } else {
                return {data.B, data.D};
            }
        }

    private:
        TransferData data;
        bool use_convert_transfer_;
        bool remove_incomplete_intersections_;
        bool handle_adaptive_refinement_;
        bool use_dual_lagrange_multiplier_;
    };
}  // namespace utopia

#endif  // UTOPIA_NEW_TRANSFER_ASSEMBLER_HPP
