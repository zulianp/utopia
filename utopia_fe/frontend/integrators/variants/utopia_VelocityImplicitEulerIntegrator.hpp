#ifndef UTOPIA_VELOCITY_IMPLICIT_EULER_INTEGRATOR_HPP
#define UTOPIA_VELOCITY_IMPLICIT_EULER_INTEGRATOR_HPP

#include "utopia_FEModelFunction.hpp"
#include "utopia_ImplicitEulerIntegrator.hpp"

#include <utility>

namespace utopia {

    template <class FunctionSpace>
    class VelocityImplicitEulerIntegrator final : public ImplicitEulerIntegrator<FunctionSpace> {
    public:
        using Super = utopia::ImplicitEulerIntegrator<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

        template <class... Args>
        VelocityImplicitEulerIntegrator(Args &&... args) : Super(std::forward<Args>(args)...) {}

        void read(Input &in) override { Super::read(in); }

        bool setup_IVP(Vector_t &x) override { return Super::setup_IVP(x); }

        bool update_IVP(const Vector_t &velocity) override {
            return Super::update_IVP(this->x_old() + this->delta_time() * velocity);
        }

        void integrate_gradient(const Vector_t &velocity, Vector_t &g) const override {
            Vector_t x = this->x_old() + this->delta_time() * velocity;
            Super::integrate_gradient(x, g);
            g *= 1. / this->delta_time();
        }

        bool time_derivative(const Vector_t &x, Vector_t &dfdt) const override {
            return Super::time_derivative(x, dfdt);
        }

        void integrate_hessian(const Vector_t &, Matrix_t &H) const override {
            const Scalar_t dt = this->delta_time();
            H *= dt;
            H += (*this->mass_matrix());
        }
    };

}  // namespace utopia

#endif  // UTOPIA_VELOCITY_IMPLICIT_EULER_INTEGRATOR_HPP
