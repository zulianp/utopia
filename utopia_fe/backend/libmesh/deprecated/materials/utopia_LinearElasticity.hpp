#ifndef UTOPIA_LINEAR_ELASTICITY_HPP
#define UTOPIA_LINEAR_ELASTICITY_HPP

#include "utopia_ElasticMaterial.hpp"
#include "utopia_fe_EDSL.hpp"

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

        void clear() override { initialized_ = false; }

        bool is_linear() const override { return true; }

        bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override {
            if (!init()) {
                return false;
            }

            hessian = constant_hessian_;
            gradient = hessian * x;
            return true;
        }

        bool assemble_gradient(const Vector &x, Vector &gradient) override {
            if (!init()) {
                return false;
            }

            gradient = constant_hessian_ * x;
            return true;
        }

        bool assemble_hessian(Matrix &hessian) override;
        bool assemble_hessian(const Vector &, Matrix &hessian) override { return assemble_hessian(hessian); }

        bool stress(const Vector &x, Vector &result) override {
            init();

            result = constant_hessian_ * x;
            result *= 1. / rescaling_;
            return true;
        }

        bool normal_stress(const UVector &x, UVector &out, const int subspace = 0) override;

        bool von_mises_stress(const UVector &x, UVector &out, const int subspace = 0) override;

        inline Scalar rescaling() const override { return rescaling_; }

        inline void rescaling(const Scalar &value) override { rescaling_ = value; }

        inline void initialize() override { init(); }

    private:
        FunctionSpaceT &V_;
        LameeParameters params_;
        bool initialized_;
        Scalar rescaling_;
        Matrix constant_hessian_;

        bool init() {
            if (initialized_) return true;

            initialized_ = assemble_hessian(constant_hessian_);
            return initialized_;
        }
    };
}  // namespace utopia

#endif  // UTOPIA_LINEAR_ELASTICITY_HPP
