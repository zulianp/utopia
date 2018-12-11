#include "utopia_QuadratureBasedAssembler.hpp"
#include "utopia_QMortarBuilder.hpp"

namespace utopia {

	QuadratureBasedAssembler::QuadratureBasedAssembler() {}
	QuadratureBasedAssembler::~QuadratureBasedAssembler() {}

	void QuadratureBasedAssembler::set_q_builder(const std::shared_ptr<QMortarBuilder> &q_builder)
	{
		this->q_builder = q_builder;
	}
}
