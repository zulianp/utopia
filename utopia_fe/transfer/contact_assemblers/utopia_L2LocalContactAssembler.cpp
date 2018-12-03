#include "utopia_L2LocalContactAssembler.hpp"

namespace utopia {

	L2LocalContactAssembler::L2LocalContactAssembler(
		const int dim,
		const double search_radius,
		const bool use_biorth)
	: q_trial_(dim), q_test_(dim), search_radius_(search_radius), use_biorth_(use_biorth)
	{}

	bool L2LocalContactAssembler::assemble(
		const Elem &trial,
		const int trial_side,
		FEType trial_type,
		const Elem &test,
		const int test_side,
		FEType test_type,
		L2LocalContactAssembler::Result &result
		)
	{
		if(!trial.has_affine_map() ||
			!test.has_affine_map()) {
			
			if(!q_builder_ || q_builder_->is_affine()) {
				q_builder_ = utopia::make_unique<WarpedContactQMortarBuilder3>(search_radius_);
			}

		} else {
			if(!q_builder_ || !q_builder_->is_affine()) {
				q_builder_ = utopia::make_unique<AffineContactQMortarBuilder3>(search_radius_);
			}
		}

		if(!q_builder_->build(
			trial,
			trial_type,
			trial_side,
			test,
			test_type,
			test_side,
			q_trial_,
			q_test_))
		{
			return false;
		}

		trial_fe_->reinit(&trial, trial_side);
		test_fe_->reinit(&test, test_side);

		const auto &gap = q_builder_->gap();
		const auto &normal = q_builder_->normal();

		if(use_biorth_) {
			// assemble_trace_biorth_weights_from_space(test_type,
			// 										 node_is_boundary_test,
			// 										 biorth_weights);

			// mortar_assemble_weighted_biorth(*trial_fe, *test_fe, biorth_weights, elemmat);
		} else {
			mortar_assemble(*trial_fe_, *test_fe_, result.coupling_matrix);
			mortar_assemble(*test_fe_,  *test_fe_, result.mass_matrix);
			integrate_scalar_function(*test_fe_, gap, result.gap);
			integrate_point_function(3, *test_fe_, normal, result.normals);

			// mortar_normal_and_gap_assemble(
			// 							   dim,
			// 							   *slave_fe,
			// 							   n_slave,
			// 							   n_master,
			// 							   plane_offset,
			// 							   surface_assemble->normals,
			// 							   surface_assemble->gap);
		}

		return false;
	}

	void L2LocalContactAssembler::print_stats(std::ostream &os) const
	{

	}

	void L2LocalContactAssembler::init_fe(
		const Elem &trial,
		FEType trial_type,
		const Elem &test,
		FEType test_type)
	{
		if(trial_fe_) return;

		trial_fe_ = libMesh::FEBase::build(trial.dim()-1, trial_type);
		test_fe_  = libMesh::FEBase::build(test.dim()-1,  test_type);
		trial_fe_->get_phi();
		test_fe_->get_phi();
		test_fe_->get_JxW();
	}

}
