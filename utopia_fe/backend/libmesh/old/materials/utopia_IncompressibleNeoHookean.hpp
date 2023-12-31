#ifndef UTOPIA_INCOMPRESSIBLE_NEOHOOKEAN_HPP
#define UTOPIA_INCOMPRESSIBLE_NEOHOOKEAN_HPP

#include "utopia.hpp"
#include "utopia_HyperElasticMaterial.hpp"
#include "utopia_fe_EDSL.hpp"
#include "utopia_libmesh_FormEval.hpp"

namespace utopia {

    template <class FunctionSpace, class Matrix, class Vector>
    class IncompressibleNeoHookean : public HyperElasticMaterial<Matrix, Vector> {
    public:
        IncompressibleNeoHookean(FunctionSpace &V, const LameeParameters &params) : V_(V), params_(params) {}

        bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override {
            // int mesh_dimension = V_.subspace(0).mesh().mesh_dimension();

            // auto mu = params_.var_mu();
            // auto lambda = params_.var_lambda();

            // auto u = trial(V_);
            // auto v = test(V_);
            // auto uk = interpolate(x, u);

            // auto F 		 = identity() + grad(uk);
            // auto F_t 	 = transpose(F);
            // auto F_inv   = inv(F);
            // auto F_inv_t = transpose(F_inv);
            // auto J       = det(F);
            // auto C       = F*F_t;
            // auto C_inv   = inv(C);

            // auto S_bar = mu * identity(mesh_dimension, mesh_dimension);

            // double d = mesh_dimension;

            // auto J_23 = power(J, -2.0/d); // to do

            // auto S_iso = S_bar + inner((-1.0 / d) * S_bar, C) * C_inv;

            // auto P1 = J_23 * (F * S_iso);

            // auto stress_lin = mu * grad(u)
            // -(lambda * logn(J) - mu) * F_inv_t * transpose(grad(u)) * F_inv_t
            // + inner(lambda * F_inv_t, grad(u)) * F_inv_t;

            // auto C_lin = (F_t * grad(u) + transpose(grad(u)) * F);

            // auto stress_lin1 = -1.0 /d * inner(C_inv, C_lin) * P1 + grad(u) * S_iso +
            //                     J_23 * F * ((-1.0 / d ) * inner(S_bar, C_lin) * C_inv +
            //                      (1.0 / d) * inner(S_bar,C) * C_inv * C_lin * C_inv);

            // auto l_form = inner(P1, grad(v)) * dX;
            // auto b_form = inner(stress_lin, grad(v)) * dX;
            // return assemble(b_form == l_form, hessian, gradient);

            assert(false);
            return false;
        }

    private:
        FunctionSpace &V_;
        LameeParameters params_;
    };
}  // namespace utopia

#endif  // UTOPIA_INC_NEOOHOOKEAN_HPP
