#ifndef UTOPIA_LINEAR_ELASTICITY_HPP
#define UTOPIA_LINEAR_ELASTICITY_HPP

#include "utopia_ElasticMaterial.hpp"
#include "ui/utopia_UIScalarSampler.hpp"
#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_fe_core.hpp"
#include "utopia_libmesh.hpp"
#include "utopia_FEFilter.hpp"
#include "utopia_FEEval_Filter.hpp"
#include "utopia_VonMisesStress.hpp"


#include <libmesh/tensor_value.h>
#include <libmesh/fe.h>

namespace utopia {

    template<class FunctionSpaceT, class Matrix, class Vector>
    class LinearElasticity final : public ElasticMaterial<Matrix, Vector> {
    public:
        LinearElasticity(FunctionSpaceT &V, const LameeParameters &params)
        : V_(V), params_(params), initialized_(false)
        {}

        bool init(Matrix &hessian)
        {
            if(initialized_) return true;
            initialized_ = assemble_hessian(hessian);
            return initialized_;
        }

        void clear() override
        {
            initialized_ = false;
        }

        bool is_linear() const override { return true; }

        
        // bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override
        bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override
        {
            if(!init(hessian)) {
                return false;
            }

            gradient = hessian * x;
            return true;
        }

        bool stress(const Vector &x, Vector &result) override {
            Matrix hessian;

            if(!assemble_hessian(hessian)) {
                return false;
            }

            result = hessian * x;
            return true;
        }

        bool normal_stress(const UVector &x, UVector &out, const int subspace = 0) override
        {
            auto u  = trial(V_);
            auto vx = test(V_[subspace]);
            
            auto mu     = params_.var_mu();
            auto lambda = params_.var_lambda();

            auto uk = interpolate(x, u);

            auto strain = transpose(grad(uk)) + grad(uk);
            auto stress = mu * strain + lambda * trace(strain) * identity();
            auto normal_stress = dot(normal(), stress * normal());
            
            UVector mass_vector;
            bool ok = utopia::assemble(surface_integral(dot(normal_stress, vx)), out); assert(ok);
            if(!ok) return false;

            utopia::assemble(surface_integral(dot(coeff(1.0), vx)), mass_vector);

            e_pseudo_inv(mass_vector, mass_vector, 1e-14);
            out = e_mul(mass_vector, out);
            return true;
        }  

        bool von_mises_stress(const UVector &x, UVector &out, const int subspace = 0) override
        {
            auto u = trial(V_);
            auto vx = test(V_[subspace]);
            
            auto mu     = params_.var_mu();
            auto lambda = params_.var_lambda();

            auto uk = interpolate(x, u);

            auto strain = transpose(grad(uk)) + grad(uk);
            auto stress = mu * strain + lambda * trace(strain) * identity();

            auto vm = filter(stress, VonMisesStress::apply);

            UVector mass_vector;
            bool ok = utopia::assemble(integral(dot(vm, vx)), out); assert(ok);
            if(!ok) return false;

            utopia::assemble(integral(dot(coeff(1.0), vx)), mass_vector);

            e_pseudo_inv(mass_vector, mass_vector, 1e-14);
            out = e_mul(mass_vector, out);
            return true;
        }

    private:
        FunctionSpaceT &V_;
        LameeParameters params_;
        bool initialized_;

        bool assemble_hessian(Matrix &hessian)
        {
            auto u = trial(V_);
            auto v = test(V_);

            auto mu     = params_.var_mu();
            auto lambda = params_.var_lambda();

            auto e_u = 0.5 * ( transpose(grad(u)) + grad(u) );
            auto e_v = 0.5 * ( transpose(grad(v)) + grad(v) );

            auto b_form = integral((2. * mu) * inner(e_u, e_v) + lambda * inner(div(u), div(v)));

            return assemble(b_form, hessian);
        }

    };
}

#endif //UTOPIA_LINEAR_ELASTICITY_HPP
