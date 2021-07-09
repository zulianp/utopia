#ifndef UTOPIA_VELOCITY_NEWMARK_INTEGRATOR_HPP
#define UTOPIA_VELOCITY_NEWMARK_INTEGRATOR_HPP

#include "utopia_FEModelFunction.hpp"
#include "utopia_NewmarkIntegrator.hpp"

#include <utility>

namespace utopia {

    template <class FunctionSpace>
    class VelocityNewmarkIntegrator final : public NewmarkIntegrator<FunctionSpace> {
    public:
        using Super = utopia::NewmarkIntegrator<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

        template <class... Args>
        VelocityNewmarkIntegrator(Args &&... args) : Super(std::forward<Args>(args)...) {}

        void read(Input &in) override { Super::read(in); }

        bool setup_IVP(Vector_t &x) override { return Super::setup_IVP(x); }

        bool update_IVP(const Vector_t &velocity) override {
            Vector_t x = this->x_old();
            update_x(velocity, x);
            velocity_old_ = velocity;
            return Super::update_IVP(x);
        }

        void integrate_gradient(const Vector_t &velocity, Vector_t &g) const override {
            Vector_t x = this->x_old();
            update_x(velocity, x);

            Super::integrate_gradient(x, g);
            g *= 2. / this->delta_time();
        }

        bool time_derivative(const Vector_t &x, Vector_t &dfdt) const override {
            return Super::time_derivative(x, dfdt);
        }

        void integrate_hessian(const Vector_t &, Matrix_t &H) const override {
            const Scalar_t dt2 = this->delta_time() * this->delta_time();
            H *= (dt2 / 4.);
            H += (*this->mass_matrix());
            this->space()->apply_constraints(H);
        }

    private:
        void update_x(const Vector_t &velocity, Vector_t &x) const {
            x += (0.5 * this->delta_time()) * (velocity_old_ + velocity);
        }

        Vector_t velocity_old_;
    };

}  // namespace utopia

#endif  // UTOPIA_VELOCITY_NEWMARK_INTEGRATOR_HPP
