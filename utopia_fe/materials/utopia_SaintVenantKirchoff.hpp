#ifndef SAINT_VENANT_KIRCHHOFF_HPP
#define SAINT_VENANT_KIRCHHOFF_HPP


#include "utopia.hpp"
#include "utopia_fe_core.hpp"
#include "utopia_HyperElasticMaterial.hpp"

namespace utopia {

	template<class FunctionSpace, class Matrix, class Vector>
	class SaintVenantKirchoff : public HyperElasticMaterial<Matrix, Vector> {
	public:

		SaintVenantKirchoff(FunctionSpace &V, const double mu, const double lambda)
		: V_(V), mu_(mu), lambda_(lambda)
		{}

		bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override
		{
			auto u = trial(V_);
			auto v = test(V_);
			auto uk = interpolate(x, u);

			auto F 		 = identity() + grad(uk);
			auto F_t 	 = transpose(F);
			auto F_inv   = inv(F);
			auto F_inv_t = transpose(F_inv);
			auto J       = det(F);
			
			auto C = F_t * F;
			auto E = 0.5 * (C - identity());
			auto S = 2.0 * mu_ * E + lambda_ * (trace(E) * identity());
			auto P = F * S;

			auto strain_lin = 0.5 * (F_t * grad(u) + transpose(grad(u)) * F);
			auto stress_lin = F * (
				2.0 * mu_ * strain_lin + lambda_ * (trace(strain_lin) * identity())
				) + grad(u) * S;


			auto l_form = inner(P, grad(v)) * dX;
			auto b_form = inner(stress_lin, grad(v)) * dX;
		
			assemble(b_form == l_form, hessian, gradient);
		}

	private:
		FunctionSpace &V_;
		double mu_, lambda_;
	};
}

#endif //SAINT_VENANT_KIRCHHOFF_HPP
