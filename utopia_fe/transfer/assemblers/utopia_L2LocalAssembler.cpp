#include "utopia_L2LocalAssembler.hpp"
#include "utopia_QMortarBuilder.hpp"

#include <cmath>
#include <queue>
#include <algorithm>
#include <sstream>
#include <numeric>

namespace utopia {
	L2LocalAssembler::L2LocalAssembler(const int dim, const bool use_biorth)
	: dim(dim),
	use_biorth(use_biorth),
	must_compute_biorth(use_biorth),
	composite_ir(dim),
	q_trial(dim),
	q_test(dim)
	{
		if(dim == 2) {
			q_builder = std::make_shared<QMortarBuilder2>();
		} else {
			assert(dim == 3);
			q_builder = std::make_shared<QMortarBuilder3>(); 
		}
	}

	L2LocalAssembler::~L2LocalAssembler()
	{}

	bool L2LocalAssembler::assemble(
		const Elem &trial,
		FEType trial_type,
		const Elem &test,
		FEType test_type,
		Matrix &mat
		)
	{
		auto trial_fe   = libMesh::FEBase::build(trial.dim(), trial_type);
		auto test_fe    = libMesh::FEBase::build(test.dim(),  test_type);
		const int order = order_for_l2_integral(dim, trial, trial_type.order, test, test_type.order);

		if(!q_builder->build(trial, trial_type, test, test_type, q_trial, q_test)) {
			return false;
		}

		init_biorth(test, test_type);
		init_fe(trial, trial_type, test, test_type);

		trial_fe->attach_quadrature_rule(&q_trial);
		trial_fe->get_phi();
		trial_fe->reinit(&trial);

		test_fe->attach_quadrature_rule(&q_test);
		test_fe->get_phi();
		test_fe->get_JxW();
		test_fe->reinit(&test);

		if(use_biorth) {
			mortar_assemble_weighted_biorth(*trial_fe, *test_fe, biorth_weights, mat);
		} else {
			mortar_assemble(*trial_fe, *test_fe, mat);
		}

		return true;
	}

	void L2LocalAssembler::init_fe(
		const Elem &trial,
		FEType trial_type,
		const Elem &test,
		FEType test_type)
	{
		if(trial_fe) return;

		trial_fe = libMesh::FEBase::build(trial.dim(), trial_type);
		test_fe  = libMesh::FEBase::build(test.dim(),  test_type);
	}

	void L2LocalAssembler::init_biorth(const Elem &test, FEType test_type)
	{
		if(!use_biorth) return;

		assemble_biorth_weights(
			test,
			test.dim(),
			test_type,
			test_type.order,
			biorth_weights);

		must_compute_biorth = false;
	}

	void L2LocalAssembler::assemble_biorth_weights(
		const libMesh::Elem &el,
		const int dim,
		const libMesh::FEType &var_type,
		const int el_order,
		libMesh::DenseMatrix<libMesh::Real> &weights)
	{
		std::unique_ptr<libMesh::FEBase> biorth_elem = libMesh::FEBase::build(dim, var_type);

		const int order = order_for_l2_integral(dim, el, el_order, el, el_order);

		libMesh::QGauss qg(dim, libMesh::Order(order));
		biorth_elem->attach_quadrature_rule(&qg);
		biorth_elem->reinit(&el);
		mortar_assemble_weights(*biorth_elem, weights);
	}
}
