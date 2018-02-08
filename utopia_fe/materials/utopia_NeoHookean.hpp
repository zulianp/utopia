#ifndef UTOPIA_NEOOHOOKEAN_HPP
#define UTOPIA_NEOOHOOKEAN_HPP

#include "utopia.hpp"
#include "utopia_fe_core.hpp"
#include "utopia_HyperElasticMaterial.hpp"

namespace utopia {

	template<class FunctionSpace, class Matrix, class Vector>
	class NeoHookean : public HyperElasticMaterial<Matrix, Vector> {
	public:

		NeoHookean(FunctionSpace &V, const double mu, const double lambda)
		: V_(V), mu_(mu), lambda_(lambda)
		{}

		bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override
		{
			Vector x_copy = x;

			auto u = trial(V_);
			auto v = test(V_);
			auto uk = interpolate(x_copy, u);

			auto F 		 = identity() + grad(uk);
			auto F_t 	 = transpose(F);
			auto F_inv   = inv(F);
			auto F_inv_t = transpose(F_inv);
			auto J       = det(F);
			
			auto P = mu_ * (F - F_inv_t) + (lambda_ * logn(J)) * F_inv_t;

			auto stress_lin = mu_ * grad(u) 
			-(lambda_ * logn(J) - mu_) * F_inv_t * transpose(grad(u)) * F_inv_t 
			+ inner(lambda_ * F_inv_t, grad(u)) * F_inv_t;

			auto l_form = inner(P, grad(v)) * dX;
			auto b_form = inner(stress_lin, grad(v)) * dX;
		
			return assemble(b_form == l_form, hessian, gradient);
		}

	private:
		FunctionSpace &V_;
		double mu_, lambda_;
	};
}

#endif //UTOPIA_NEOOHOOKEAN_HPP

