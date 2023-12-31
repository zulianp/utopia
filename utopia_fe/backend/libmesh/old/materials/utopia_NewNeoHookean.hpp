#ifndef UTOPIA_NEW_NEO_HOOKEAN_HPP
#define UTOPIA_NEW_NEO_HOOKEAN_HPP

#include "utopia.hpp"
#include "utopia_HyperElasticMaterial.hpp"
#include "utopia_fe_EDSL.hpp"

namespace utopia {

    template <class FunctionSpace, class Matrix, class Vector>
    class NewNeoHookean : public HyperElasticMaterial<Matrix, Vector> {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);

        NewNeoHookean(FunctionSpace &V, const LameeParameters &params) : V_(V), params_(params), rescaling_(1.0) {}

        bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override;
        bool assemble_hessian(const Vector &x, Matrix &hessian) override;
        bool assemble_gradient(const Vector &x, Vector &gradient) override;

        bool stress(const Vector &x, Vector &result) override;

        inline Scalar rescaling() const override { return rescaling_; }

        inline void rescaling(const Scalar &value) override { rescaling_ = value; }

    private:
        FunctionSpace &V_;
        LameeParameters params_;
        Scalar rescaling_;
    };
}  // namespace utopia

#endif  // UTOPIA_NEOOHOOKEAN_HPP
