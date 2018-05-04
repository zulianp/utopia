#ifndef UTOPIA_STABILIZED_MATERIAL_HPP
#define UTOPIA_STABILIZED_MATERIAL_HPP

#include "utopia.hpp"
#include "utopia_fe_core.hpp"
#include "utopia_ElasticMaterial.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"

namespace utopia {

	template<class FunctionSpace, class Matrix, class Vector>
	class StabilizedMaterial : public ElasticMaterial<Matrix, Vector> {
	public:

		StabilizedMaterial(
			FunctionSpace &V,
			const double stabilization_mag,
			const std::shared_ptr<ElasticMaterial<Matrix, Vector>> &material)
		: V_(V), stabilization_mag_(stabilization_mag), material_(material)
		{}

		bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override
		{
			if(!material_->assemble_hessian_and_gradient(x, hessian_, gradient)) {
				std::cerr << "decorated material failed to assemble" << std::endl;
				return false;
			}

			if(empty(lapl_)) {
				auto u = trial(V_);
				auto v = test(V_);

				auto b_form = inner(u, v) * dX;
				utopia::assemble(b_form, lapl_);
				lapl_ *= stabilization_mag_;
			}

			hessian = hessian_ + lapl_;
			return true;
		}

	private:
		FunctionSpace &V_;
		double stabilization_mag_;
		std::shared_ptr<ElasticMaterial<Matrix, Vector>> material_;

		Matrix hessian_;
		Matrix lapl_;
	};
}

#endif //UTOPIA_NEOOHOOKEAN_HPP

