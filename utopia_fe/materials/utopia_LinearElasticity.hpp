#ifndef UTOPIA_LINEAR_ELASTICITY_HPP
#define UTOPIA_LINEAR_ELASTICITY_HPP

#include "utopia_ElasticMaterial.hpp"
#include "utopia_fe_core.hpp"

#include "utopia_Integrators.hpp"
#include "utopia_LameeParameters.hpp"
#include "utopia_libmesh_FormEval.hpp"

namespace utopia {

    template <class FunctionSpaceT, class Matrix, class Vector>
    class LinearElasticity final
        : public ElasticMaterial<Matrix, Vector>  //, public ModularEquationIntegrator<FunctionSpaceT>
    {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);

        LinearElasticity(FunctionSpaceT &V, const LameeParameters &params);
        ~LinearElasticity();

        bool init(Matrix &hessian) {
            if (initialized_) return true;

            initialized_ = assemble_hessian(hessian);
            return initialized_;
        }

        void clear() override { initialized_ = false; }

        bool is_linear() const override { return true; }

        bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override {
            if (!init(hessian)) {
                return false;
            }

            gradient = hessian * x;
            return true;
        }

        bool stress(const Vector &x, Vector &result) override {
            Matrix hessian;

            if (!assemble_hessian(hessian)) {
                return false;
            }

            result = hessian * x;
            result *= 1. / rescaling_;
            return true;
        }

        bool normal_stress(const UVector &x, UVector &out, const int subspace = 0) override;

        bool von_mises_stress(const UVector &x, UVector &out, const int subspace = 0) override;

        inline Scalar rescaling() const override { return rescaling_; }

        inline void rescaling(const Scalar &value) override { rescaling_ = value; }

    private:
        FunctionSpaceT &V_;
        LameeParameters params_;
        bool initialized_;
        Scalar rescaling_;

        // inline void init_integrators(const UVector &x)
        // {
        //     init_bilinear_integrator();
        //     init_linear_integrator(x);
        // }

        // void init_bilinear_integrator();

        // void init_linear_integrator(const UVector &x);

        bool assemble_hessian(Matrix &hessian);
    };
}  // namespace utopia

#endif  // UTOPIA_LINEAR_ELASTICITY_HPP
