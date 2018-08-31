#ifndef UTOPIA_POUROUS_MEDIA_TO_FRACTURE_TRANSFER_HPP
#define UTOPIA_POUROUS_MEDIA_TO_FRACTURE_TRANSFER_HPP

#include "utopia_libmesh.hpp"
#include "utopia_TransferAssembler.hpp"

#include <memory>

namespace utopia {
	class MeshTransferOperator final : public TransferOperator {
	public:
		using SparseMatrix  = utopia::DSMatrixd;
		using Vector 		= utopia::DVectord;
		using MeshBase      = libMesh::MeshBase;
		using DofMap        = libMesh::DofMap;

		MeshTransferOperator(
			const std::shared_ptr<MeshBase> &from_mesh,
			const std::shared_ptr<DofMap>   &from_dofs,
			const std::shared_ptr<MeshBase> &to_mesh,
			const std::shared_ptr<DofMap>   &to_dofs,
			const TransferOptions &opts = TransferOptions()
		) : 
		from_mesh(from_mesh),
		from_dofs(from_dofs),
		to_mesh(to_mesh),
		to_dofs(to_dofs),
		opts(opts)
		{}

		//@brief operator_type \in \{ INTERPOLATION| L2_PROJECTION| PSEUDO_L2_PROJECTION | APPROX_L2_PROJECTION \}
		bool initialize(const TransferOperatorType operator_type = utopia::INTERPOLATION);
		bool initialize(const std::string operator_type);

		inline void apply(const Vector &from, Vector &to) const
		{
			operator_->apply(from, to);
		}

		inline void apply_transpose(const Vector &from, Vector &to) const
		{
			operator_->apply_transpose(from, to);
		}

		inline void describe(std::ostream &os) const
		{
			if(operator_) {
				operator_->describe(os);
			}
		}

	private:
		std::shared_ptr<MeshBase> from_mesh;
		std::shared_ptr<DofMap>   from_dofs;
		std::shared_ptr<MeshBase> to_mesh;
		std::shared_ptr<DofMap>   to_dofs;
		TransferOptions opts;

		std::shared_ptr<TransferOperator> operator_;

	};

	using PourousMediaToFractureTransfer = MeshTransferOperator;
}

#endif //UTOPIA_POUROUS_MEDIA_TO_FRACTURE_TRANSFER_HPP
