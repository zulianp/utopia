#ifndef UTOPIA_POUROUS_MEDIA_TO_FRACTURE_TRANSFER_R_HPP
#define UTOPIA_POUROUS_MEDIA_TO_FRACTURE_TRANSFER_R_HPP

#include "utopia_libmesh.hpp"
#include "utopia_TransferAssemblerR.hpp"
#include "utopia_Path.hpp"


#include <memory>

namespace utopia {
	class MeshTransferOperatorBidirectional final : public TransferOperatorR {
	public:
		using SparseMatrix  = utopia::USparseMatrix;
		using Vector 		= utopia::UVector;
		using MeshBase      = libMesh::MeshBase;
		using DofMap        = libMesh::DofMap;

		MeshTransferOperatorBidirectional(
			const std::shared_ptr<MeshBase> &from_mesh,
			const std::shared_ptr<DofMap>   &from_dofs,
			const std::shared_ptr<DofMap>   &from_dofs_r,
			const std::shared_ptr<MeshBase> &to_mesh,
			const std::shared_ptr<DofMap>   &to_dofs,
			const std::shared_ptr<DofMap>   &to_dofs_r,
			const TransferOptionsR &opts = TransferOptionsR()
		) : 
		from_mesh(from_mesh),
		from_dofs(from_dofs),
	    from_dofs_r(from_dofs_r),
		to_mesh(to_mesh),
		to_dofs(to_dofs),
		to_dofs_r(to_dofs_r),
		opts(opts)
		{}

		//@brief operator_type \in \{ INTERPOLATION| L2_PROJECTION| PSEUDO_L2_PROJECTION | APPROX_L2_PROJECTION \}
		bool initialize(const TransferOperatorRType operator_type = utopia::BIDIRECTIONAL_L2_PROJECTION);
		bool initialize(const std::string operator_type);

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

		void set_normalize_rows(const bool val)
		{

		}

	private:
		std::shared_ptr<MeshBase> from_mesh;
		std::shared_ptr<DofMap>   from_dofs;
		std::shared_ptr<DofMap>   from_dofs_r;
		std::shared_ptr<MeshBase> to_mesh;
		std::shared_ptr<DofMap>   to_dofs;
		std::shared_ptr<DofMap>   to_dofs_r;
		TransferOptionsR opts;

		std::shared_ptr<TransferOperatorR> operator_;
		bool normalize_rows_;

	};

	using PourousMediaToFractureTransfer = MeshTransferOperatorBidirectional;
}

#endif //UTOPIA_POUROUS_MEDIA_TO_FRACTURE_TRANSFER_HPP
