#include "utopia_MeshTransferOperatorBidirectional.hpp"

#include "utopia_L2LocalAssembler.hpp"
#include "utopia_ApproxL2LocalAssembler.hpp"
#include "utopia_InterpolationLocalAssembler.hpp"
#include "utopia_Local2Global.hpp"
#include "utopia_BidirectionalL2LocalAssembler.hpp"

#include <map>

namespace utopia {

	static const std::map<std::string, TransferOperatorRType> &get_str_to_type()
	{
		static std::map<std::string, TransferOperatorRType> types;

		if(types.empty()) {
			types["BIDIRECTIONAL_L2_PROJECTION"] = BIDIRECTIONAL_L2_PROJECTION;
			types["BIDIRECTIONAL_PSEUDO_L2_PROJECTION"] = BIDIRECTIONAL_PSEUDO_L2_PROJECTION;

			//other way of writing them
			types["bidirectional-l2-projection"] = BIDIRECTIONAL_L2_PROJECTION;
			types["bidirectional-pseudo-l2-projection"] = BIDIRECTIONAL_PSEUDO_L2_PROJECTION;
		}

		return types;
	}

	bool MeshTransferOperatorBidirectional::initialize(const std::string operator_type)
	{
		const auto &m = get_str_to_type(); 

		auto it = m.find(operator_type);

		if(it == m.end()) {
			return initialize(BIDIRECTIONAL_PSEUDO_L2_PROJECTION);
		} else {
			return initialize(it->second);
		}
	}

	bool MeshTransferOperatorBidirectional::initialize(const TransferOperatorRType operator_type)
	{
		std::shared_ptr<LocalAssembler> assembler;

		bool use_interpolation = false;
		bool use_biorth        = false;
		bool is_bidirectonal   = false;

		switch(operator_type) {
			

			case BIDIRECTIONAL_L2_PROJECTION:
			{
				std::cout << "[Status] using bi l2 projection" << std::endl;
				assembler = std::make_shared<BidirectionalL2LocalAssembler>(from_mesh->mesh_dimension(), false, true);
				is_bidirectonal = true;
				break;
			}

			case BIDIRECTIONAL_PSEUDO_L2_PROJECTION:
			{
				std::cout << "[Status] using bi pseudo l2 projection" << std::endl;
				assembler = std::make_shared<BidirectionalL2LocalAssembler>(from_mesh->mesh_dimension(), true, false);
				use_biorth = true;
				is_bidirectonal = true;
				break;
			}

			default:
			{
				assert(false);
				return false;
			}
		}

		auto local2global = std::make_shared<Local2Global>(use_interpolation);
		TransferAssemblerR transfer_assembler(assembler, local2global);

		std::vector< std::shared_ptr<SparseMatrix> > mats;
		if(!transfer_assembler.assemble(from_mesh, from_dofs, from_dofs_r, to_mesh, to_dofs, to_dofs_r, mats, opts)) {
			return false;
		}

		if(is_bidirectonal) {
			if(use_biorth) {
				auto forward = std::make_shared<PseudoL2TransferOperatorR>();
				forward->init_from_coupling_operator(*mats[0]);

				auto backward = std::make_shared<PseudoL2TransferOperatorR>();
				backward->init_from_coupling_operator(*mats[1]);
				operator_ = std::make_shared<BidirectionalOperator>(forward, backward);

			} else {
				auto forward = std::make_shared<L2TransferOperatorR>(mats[0], mats[1], std::make_shared<Factorization<USparseMatrix, UVector>>());
				forward->fix_mass_matrix_operator();

				auto backward = std::make_shared<L2TransferOperatorR>(mats[2], mats[3], std::make_shared<Factorization<USparseMatrix, UVector>>());
				backward->fix_mass_matrix_operator();
				operator_ = std::make_shared<BidirectionalOperator>(forward, backward);
			}

		} 

		operator_->describe(std::cout);
		return true;
	}
}
