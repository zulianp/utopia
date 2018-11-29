#ifndef UTOPIA_POUROUS_MEDIA_TO_FRACTURE_TRANSFER_HPP
#define UTOPIA_POUROUS_MEDIA_TO_FRACTURE_TRANSFER_HPP

#include "utopia_libmesh.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_Path.hpp"


#include <memory>

namespace utopia {
	class MeshTransferOperator final : public TransferOperator {
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
		) : 
		from_mesh(from_mesh),
		from_dofs(from_dofs),
		to_mesh(to_mesh),
		to_dofs(to_dofs),
		opts(opts),
		normalize_rows_(true),
		tol_(1e-14)
		{}

		//@brief operator_type \in \{ INTERPOLATION| L2_PROJECTION| PSEUDO_L2_PROJECTION | APPROX_L2_PROJECTION \}
		bool initialize(const TransferOperatorType operator_type = utopia::INTERPOLATION);
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


		void set_tol(const double val)
		{
			tol_ = val;
		}

	private:
		std::shared_ptr<MeshBase> from_mesh;
		std::shared_ptr<DofMap>   from_dofs;
		std::shared_ptr<MeshBase> to_mesh;
		std::shared_ptr<DofMap>   to_dofs;
		TransferOptions opts;

		std::shared_ptr<TransferOperator> operator_;
		bool normalize_rows_;
		double tol_;

	};

	using PourousMediaToFractureTransfer = MeshTransferOperator;
}

#endif //UTOPIA_POUROUS_MEDIA_TO_FRACTURE_TRANSFER_HPP
