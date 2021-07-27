#ifndef UTOPIA_VELOCITY_IMPLICIT_EULER_INTEGRATOR_HPP
#define UTOPIA_VELOCITY_IMPLICIT_EULER_INTEGRATOR_HPP

#include "utopia_FEModelFunction.hpp"
#include "utopia_ImplicitEulerIntegrator.hpp"
#include "utopia_VelocityVariant.hpp"

#include <utility>

namespace utopia {

    template <class FunctionSpace>
    class VelocityImplicitEulerIntegrator final : public ImplicitEulerIntegrator<FunctionSpace>,
                                                  public VelocityVariant<FunctionSpace> {
    public:
        using Super = utopia::ImplicitEulerIntegrator<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

        template <class... Args>
        VelocityImplicitEulerIntegrator(Args &&... args) : Super(std::forward<Args>(args)...) {}

        const Vector_t &solution() const override { return this->x_old(); }

        void read(Input &in) override { Super::read(in); }

        bool setup_IVP(Vector_t &x) override { return Super::setup_IVP(x); }

        bool update_IVP(const Vector_t &velocity) override {
            Vector_t x;
            update_x(velocity, x);
            return Super::update_IVP(x);
        }

        bool gradient(const Vector_t &velocity, Vector_t &g) const override {
            Vector_t x;
            update_x(velocity, x);
            return Super::gradient(x, g);
        }

        bool hessian(const Vector_t &velocity, Matrix_t &H) const override {
            Vector_t x;
            update_x(velocity, x);
            return Super::hessian(x, H);
        }

        void integrate_gradient(const Vector_t &x, Vector_t &g) const override {
            Super::integrate_gradient(x, g);
            g *= 1. / this->delta_time();
        }

        bool time_derivative(const Vector_t &x, Vector_t &dfdt) const override {
            return Super::time_derivative(x, dfdt);
        }

        void integrate_hessian(const Vector_t &x, Matrix_t &H) const override { Super::integrate_hessian(x, H); }

    private:
        void update_x(const Vector_t &velocity, Vector_t &x) const {
            x = this->x_old();
            x += this->delta_time() * velocity;
            this->space()->apply_constraints(x);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_VELOCITY_IMPLICIT_EULER_INTEGRATOR_HPP
