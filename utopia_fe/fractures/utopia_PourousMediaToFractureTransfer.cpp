#include "utopia_PourousMediaToFractureTransfer.hpp"

#include "utopia_L2LocalAssembler.hpp"
#include "utopia_InterpolationLocalAssembler.hpp"
#include "utopia_Local2Global.hpp"


namespace utopia {

	bool PourousMediaToFractureTransfer::initialize(const TransferOperatorType operator_type)
	{
		std::shared_ptr<LocalAssembler> assembler;

		bool use_interpolation = operator_type == INTERPOLATION;
		bool use_biorth = operator_type == PSEUDO_L2_PROJECTION;

		if(use_interpolation) {
			std::cout << "[Status] using interpolation" << std::endl;
			assembler = std::make_shared<InterpolationLocalAssembler>(from_mesh->mesh_dimension());
		} else {
			std::cout << "[Status] using projection" << std::endl;
			assembler = std::make_shared<L2LocalAssembler>(from_mesh->mesh_dimension(), use_biorth);
		}

		auto local2global = std::make_shared<Local2Global>(use_interpolation);
		TransferAssembler transfer_assembler(assembler, local2global);

		auto B = std::make_shared<SparseMatrix>();

		if(!transfer_assembler.assemble(from_mesh, from_dofs, to_mesh, to_dofs, *B, opts)) {
			return false;
		}

		if(use_interpolation) {
			auto interpolation_operator = std::make_shared<Interpolator>(B);
			//only necessary for non-conforming mesh in parallel
			interpolation_operator->normalize_rows();
			operator_ = interpolation_operator;
		} else if(use_biorth) {
			auto pseudo_l2_operator = std::make_shared<PseudoL2TransferOperator>();
			pseudo_l2_operator->init_from_coupling_operator(*B);
			operator_ = pseudo_l2_operator;
		} else {
			assert(false && "must assemble mass matrix. Not implemented yet");
			auto D = std::make_shared<SparseMatrix>();
			auto l2_operator = std::make_shared<L2TransferOperator>(B, D);
			operator_ = l2_operator;
		}

		operator_->describe(std::cout);
		return true;
	}
}
