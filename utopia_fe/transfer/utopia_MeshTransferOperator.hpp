#ifndef UTOPIA_MESH_TRANSFER_OPERATOR_HPP
#define UTOPIA_MESH_TRANSFER_OPERATOR_HPP

#include "utopia_libmesh.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_Path.hpp"
#include "utopia_Input.hpp"


#include <memory>
#include <functional>

namespace utopia {
    class MeshTransferOperator final : public TransferOperator, public Configurable {
    public:
        using SparseMatrix  = utopia::USparseMatrix;
        using Vector 		= utopia::UVector;
        using MeshBase      = libMesh::MeshBase;
        using DofMap        = libMesh::DofMap;

        MeshTransferOperator(
            const std::shared_ptr<MeshBase> &from_mesh,
            const std::shared_ptr<DofMap>   &from_dofs,
            const std::shared_ptr<MeshBase> &to_mesh,
            const std::shared_ptr<DofMap>   &to_dofs,
            const TransferOptions &opts = TransferOptions()
        );

        ~MeshTransferOperator();

        void read(Input &is) override;

        //@brief operator_type \in \{ INTERPOLATION| L2_PROJECTION| PSEUDO_L2_PROJECTION | APPROX_L2_PROJECTION \}
        bool initialize(const TransferOperatorType operator_type = utopia::INTERPOLATION);
        bool initialize(const std::string operator_type);
        bool assemble();

        inline void apply(const Vector &from, Vector &to) const override
        {
            assert(operator_);
            operator_->apply(from, to);
        }

        inline void apply_transpose(const Vector &from, Vector &to) const override
        {
            assert(operator_);
            operator_->apply_transpose(from, to);
        }

        inline void describe(std::ostream &os) const override
        {
            if(operator_) {
                operator_->describe(os);
            }
        }

        bool write(const Path &path) const override
        {
            if(operator_) {
                return operator_->write(path);
            }

            return false;
        }

        template<class AlgebraicOperator>
        inline std::shared_ptr<AlgebraicOperator> get() const
        {
            return std::dynamic_pointer_cast<AlgebraicOperator>(operator_);
        }

        void set_tol(const double val);

    private:
        std::shared_ptr<MeshBase> from_mesh;
        std::shared_ptr<MeshBase> filtered_from_mesh;
        std::shared_ptr<DofMap>   from_dofs;
        std::shared_ptr<MeshBase> to_mesh;
        std::shared_ptr<MeshBase> filtered_to_mesh;
        std::shared_ptr<DofMap>   to_dofs;
        TransferOptions opts;

        std::shared_ptr<TransferOperator> operator_;

        class Params;
        std::unique_ptr<Params> params_;

        //strategies
        bool set_up_l2_projection();
        bool set_up_pseudo_l2_projection();
        bool set_up_interpolation();
        bool set_up_approx_l2_projection();
        bool set_up_bidirectional_transfer();
        bool set_up_bidirectional_pseudo_transfer();
        static std::unique_ptr<LinearSolver<USparseMatrix, UVector> > new_solver();

        std::map<TransferOperatorType, std::function<bool()>> assembly_strategies_;

        /////////////
        std::shared_ptr<libMesh::MeshBase> get_filtered_from_mesh();
        std::shared_ptr<libMesh::MeshBase> get_filtered_to_mesh();


        void assemble_mass_matrix(
            const libMesh::MeshBase &mesh,
            const libMesh::DofMap &dof_map,
            const int var,
            const int n_tensor,
            USparseMatrix &mat
            ) const;
    };

}

#endif //UTOPIA_MESH_TRANSFER_OPERATOR_HPP
