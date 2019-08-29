#ifndef UTOPIA_NEOHOOKEAN_HPP
#define UTOPIA_NEOHOOKEAN_HPP

#include "utopia.hpp"
#include "utopia_fe_core.hpp"
#include "utopia_HyperElasticMaterial.hpp"

namespace utopia {

    template<class FunctionSpace, class Matrix, class Vector>
    class NeoHookean : public HyperElasticMaterial<Matrix, Vector> {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);

        NeoHookean(FunctionSpace &V, const LameeParameters &params)
        : V_(V), params_(params), rescaling_(1.0)
        {}

        bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override
        {
            auto mu     = params_.var_mu();
            auto lambda = params_.var_lambda();

            auto u  = trial(V_);
            auto v  = test(V_);
            auto uk = interpolate(x, u);

            auto F 		 = identity() + grad(uk);
            auto F_t 	 = transpose(F);
            auto F_inv   = inv(F);
            auto F_inv_t = transpose(F_inv);
            auto J       = det(F);

            auto P = (rescaling_ * mu) * (F - F_inv_t) + ((rescaling_ * lambda) * logn(J)) * F_inv_t;

            auto stress_lin = (rescaling_ * mu) * grad(u)
            - (rescaling_ * (lambda * logn(J) - mu)) * F_inv_t * transpose(grad(u)) * F_inv_t
            + inner(lambda * F_inv_t, grad(u)) * F_inv_t;

            auto l_form = inner(P, grad(v)) * dX;
            auto b_form = inner(stress_lin, grad(v)) * dX;

            return assemble(b_form == l_form, hessian, gradient);
        }
        
        bool stress(const Vector &x, Vector &result) override {

            auto mu     = params_.var_mu();
            auto lambda = params_.var_lambda();

            auto u  = trial(V_);
            auto v  = test(V_);
            auto uk = interpolate(x, u);

            auto F 		 = identity() + grad(uk);
            auto F_inv   = inv(F);
            auto F_inv_t = transpose(F_inv);
            auto J       = det(F);

            auto P = (rescaling_ * mu) * (F - F_inv_t) + ((rescaling_ * lambda) * logn(J)) * F_inv_t;

            auto l_form = inner(P, grad(v)) * dX;
            return assemble(l_form, result);
        }

        inline Scalar rescaling() const override
        {
            return rescaling_;
        }

        inline void rescaling(const Scalar &value) override {
            rescaling_ = value;
        }


    private:
        FunctionSpace &V_;
        LameeParameters params_;
        Scalar rescaling_;
    };
}

#endif //UTOPIA_NEOOHOOKEAN_HPP

