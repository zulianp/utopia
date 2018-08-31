#include "utopia_PourousMediaToFractureTransfer.hpp"

#include "utopia_L2LocalAssembler.hpp"
#include "utopia_ApproxL2LocalAssembler.hpp"
#include "utopia_InterpolationLocalAssembler.hpp"
#include "utopia_Local2Global.hpp"


namespace utopia {

	static const std::map<std::string, TransferOperatorType> &get_str_to_type()
	{
		static std::map<std::string, TransferOperatorType> types;

		if(types.empty()) {
			types["INTERPOLATION"] 		  = INTERPOLATION;
			types["L2_PROJECTION"] 		  = L2_PROJECTION;
			types["PSEUDO_L2_PROJECTION"] = PSEUDO_L2_PROJECTION;
			types["APPROX_L2_PROJECTION"] = APPROX_L2_PROJECTION;
		}

		return types;
	}

	bool MeshTransferOperator::initialize(const std::string operator_type)
	{
		const auto &m = get_str_to_type(); 

		auto it = m.find(operator_type);

		if(it == m.end()) {
			return initialize(PSEUDO_L2_PROJECTION);
		} else {
			return initialize(it->second);
		}
	}

	bool MeshTransferOperator::initialize(const TransferOperatorType operator_type)
	{
		std::shared_ptr<LocalAssembler> assembler;

		bool use_interpolation = false;
		bool use_biorth        = false;

		switch(operator_type) {
			case INTERPOLATION:
			{
				std::cout << "[Status] using interpolation" << std::endl;
				assembler = std::make_shared<InterpolationLocalAssembler>(from_mesh->mesh_dimension());
				use_interpolation = true;
				break;
			}

			case L2_PROJECTION:
			{	
				std::cout << "[Status] using l2 projection" << std::endl;
				assembler = std::make_shared<L2LocalAssembler>(from_mesh->mesh_dimension(), false, true);
				break;
			}

			case PSEUDO_L2_PROJECTION:
			{
				std::cout << "[Status] using pseudo l2 projection" << std::endl;
				assembler = std::make_shared<L2LocalAssembler>(from_mesh->mesh_dimension(), true, false);
				use_biorth = true;
				break;
			}

			case APPROX_L2_PROJECTION:
			{
				std::cout << "[Status] using approx l2 projection" << std::endl;
				assembler = std::make_shared<ApproxL2LocalAssembler>(from_mesh->mesh_dimension());
				break;
			}

			default:
			{
				assert(false);
				return false;
			}
		}

		auto local2global = std::make_shared<Local2Global>(use_interpolation);
		TransferAssembler transfer_assembler(assembler, local2global);

		std::vector< std::shared_ptr<SparseMatrix> > mats;
		if(!transfer_assembler.assemble(from_mesh, from_dofs, to_mesh, to_dofs, mats, opts)) {
			return false;
		}

		if(use_interpolation) {
			auto interpolation_operator = std::make_shared<Interpolator>(mats[0]);
			//only necessary for non-conforming mesh in parallel
			interpolation_operator->normalize_rows();
			operator_ = interpolation_operator;
		} else if(use_biorth) {
			auto pseudo_l2_operator = std::make_shared<PseudoL2TransferOperator>();
			pseudo_l2_operator->init_from_coupling_operator(*mats[0]);
			operator_ = pseudo_l2_operator;
		} else {
			auto l2_operator = std::make_shared<L2TransferOperator>(mats[0], mats[1], std::make_shared<Factorization<DSMatrixd, DVectord>>());
			l2_operator->fix_mass_matrix_operator();
			operator_ = l2_operator;
		}

		operator_->describe(std::cout);
		return true;
	}
}
