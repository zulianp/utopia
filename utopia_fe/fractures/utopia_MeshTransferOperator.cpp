#include "utopia_MeshTransferOperator.hpp"

#include "utopia_L2LocalAssembler.hpp"
#include "utopia_ApproxL2LocalAssembler.hpp"
#include "utopia_InterpolationLocalAssembler.hpp"
#include "utopia_Local2Global.hpp"
#include "utopia_BidirectionalL2LocalAssembler.hpp"

#include <map>

namespace utopia {

	static const std::map<std::string, TransferOperatorType> &get_str_to_type()
	{
		static std::map<std::string, TransferOperatorType> types;

		if(types.empty()) {
			types["INTERPOLATION"] 		  = INTERPOLATION;
			types["L2_PROJECTION"] 		  = L2_PROJECTION;
			types["PSEUDO_L2_PROJECTION"] = PSEUDO_L2_PROJECTION;
			types["APPROX_L2_PROJECTION"] = APPROX_L2_PROJECTION;
			types["BIDIRECTIONAL_L2_PROJECTION"] = BIDIRECTIONAL_L2_PROJECTION;
			types["BIDIRECTIONAL_PSEUDO_L2_PROJECTION"] = BIDIRECTIONAL_PSEUDO_L2_PROJECTION;

			//other way of writing them
			types["interpolation"] 		  = INTERPOLATION;
			types["l2-projection"] 		  = L2_PROJECTION;
			types["pseudo-l2-projection"] = PSEUDO_L2_PROJECTION;
			types["approx-l2-projection"] = APPROX_L2_PROJECTION;
			types["bidirectional-l2-projection"] = BIDIRECTIONAL_L2_PROJECTION;
			types["bidirectional-pseudo-l2-projection"] = BIDIRECTIONAL_PSEUDO_L2_PROJECTION;
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

	inline static void assemble_mass_matrix(
		const libMesh::MeshBase &mesh,
		const libMesh::DofMap &dof_map, 
		const int var,
		const int n_tensor,
		USparseMatrix &mat
		)
	{	
		auto dim = mesh.mesh_dimension();
		auto fe_type = dof_map.variable_type(var);
		std::vector<libMesh::dof_id_type> indices;

		SizeType nnz_x_row = 0;

		if(!dof_map.get_n_nz().empty()) {
			nnz_x_row = 
			*std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()) + 
			*std::max_element(dof_map.get_n_oz().begin(), dof_map.get_n_oz().end());
		}

		mat = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);

		{
			Write<USparseMatrix> w_m(mat, utopia::GLOBAL_ADD);
			libMesh::DenseMatrix<libMesh::Real> el_mat;

			auto fe = libMesh::FEBase::build(dim, fe_type);
			libMesh::QGauss qrule(dim, fe_type.default_quadrature_order());
			fe->attach_quadrature_rule(&qrule);
			
			auto &phi = fe->get_phi();
			auto &JxW = fe->get_JxW();
			
			for(auto it = elements_begin(mesh); it != elements_end(mesh); ++it) {
				fe->reinit(*it);
				dof_map.dof_indices(*it, indices);

				auto n_shape_functions = phi.size();
				el_mat.resize(n_shape_functions, n_shape_functions);
				el_mat.zero();

				for(unsigned int qp = 0; qp < qrule.n_points(); qp++) {
					for(unsigned int i = 0; i < n_shape_functions; i++) {
						for(unsigned int j = 0; j < n_shape_functions; j++) {
							auto value = JxW[qp]*phi[i][qp]*phi[j][qp];

							for(unsigned int k = 0; k < dim; k++) {
								el_mat(i + k * phi.size(), j + k * phi.size()) += value;
							}
						}
					}
				}

				add_matrix(el_mat, indices, indices, mat);
			}
		}

	}

	bool MeshTransferOperator::initialize(const TransferOperatorType operator_type)
	{
		std::shared_ptr<LocalAssembler> assembler;

		bool use_interpolation = false;
		bool use_biorth        = false;
		bool is_bidirectonal   = false;

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

			case BIDIRECTIONAL_L2_PROJECTION:
			{
				std::cout << "[Status] using bi l2 projection" << std::endl;
				assembler = std::make_shared<BidirectionalL2LocalAssembler>(from_mesh->mesh_dimension(), false, !bi_operator_mass_mat_outside_);
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
		TransferAssembler transfer_assembler(assembler, local2global);

		std::vector< std::shared_ptr<SparseMatrix> > mats;
		if(!transfer_assembler.assemble(from_mesh, from_dofs, to_mesh, to_dofs, mats, opts)) {
			return false;
		}

		if(is_bidirectonal) {

			if(bi_operator_mass_mat_outside_) {
				auto mass_mat_from = std::make_shared<USparseMatrix>();
				auto mass_mat_to   = std::make_shared<USparseMatrix>();

				assemble_mass_matrix(*from_mesh, *from_dofs, opts.from_var_num, opts.n_var, *mass_mat_from);
				assemble_mass_matrix(*to_mesh,   *to_dofs,   opts.to_var_num,   opts.n_var, *mass_mat_to);

				if(use_biorth) {
					auto forward = std::make_shared<PseudoL2TransferOperator>();
					forward->init_from_coupling_and_mass_operator(*mats[0], *mass_mat_to);

					auto backward = std::make_shared<PseudoL2TransferOperator>();
					backward->init_from_coupling_and_mass_operator(*mats[1], *mass_mat_from);
					operator_ = std::make_shared<BidirectionalOperator>(forward, backward);

				} else {
					auto forward = std::make_shared<L2TransferOperator>(mats[0], mass_mat_to, true, std::make_shared<Factorization<USparseMatrix, UVector>>());
					forward->fix_mass_matrix_operator();

					auto backward = std::make_shared<L2TransferOperator>(mats[1], mass_mat_from, true, std::make_shared<Factorization<USparseMatrix, UVector>>());
					backward->fix_mass_matrix_operator();
					operator_ = std::make_shared<BidirectionalOperator>(forward, backward);
				}

			} else {
				if(use_biorth) {
					auto forward = std::make_shared<PseudoL2TransferOperator>();
					forward->init_from_coupling_operator(*mats[0]);

					auto backward = std::make_shared<PseudoL2TransferOperator>();
					backward->init_from_coupling_operator(*mats[1]);
					operator_ = std::make_shared<BidirectionalOperator>(forward, backward);

				} else {
					auto forward = std::make_shared<L2TransferOperator>(mats[0], mats[1], std::make_shared<Factorization<USparseMatrix, UVector>>());
					forward->fix_mass_matrix_operator();

					auto backward = std::make_shared<L2TransferOperator>(mats[2], mats[3], std::make_shared<Factorization<USparseMatrix, UVector>>());
					backward->fix_mass_matrix_operator();
					operator_ = std::make_shared<BidirectionalOperator>(forward, backward);
				}
			}

		} else {
			if(use_interpolation) {
				auto interpolation_operator = std::make_shared<Interpolator>(mats[0]);
				//only necessary for non-conforming mesh in parallel
				interpolation_operator->normalize_rows();
				operator_ = interpolation_operator;
			} else if(use_biorth) {

				if(normalize_rows_) {
					auto pseudo_l2_operator = std::make_shared<PseudoL2TransferOperator>();
					pseudo_l2_operator->init_from_coupling_operator(*mats[0]);
					operator_ = pseudo_l2_operator;
				} else {
					operator_ = std::make_shared<PseudoL2TransferOperator>(mats[0]);
				}

			} else {
				auto l2_operator = std::make_shared<L2TransferOperator>(mats[0], mats[1], std::make_shared<Factorization<USparseMatrix, UVector>>());
				l2_operator->fix_mass_matrix_operator(tol_);
				operator_ = l2_operator;
			}
		}

		operator_->describe(std::cout);
		return true;
	}
}
