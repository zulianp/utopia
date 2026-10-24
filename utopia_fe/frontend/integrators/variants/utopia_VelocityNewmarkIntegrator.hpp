#ifndef UTOPIA_VELOCITY_NEWMARK_INTEGRATOR_HPP
#define UTOPIA_VELOCITY_NEWMARK_INTEGRATOR_HPP

#include "utopia_Tracer.hpp"

#include "utopia_FEModelFunction.hpp"
#include "utopia_NewmarkIntegrator.hpp"
#include "utopia_VelocityVariant.hpp"

#include <utility>

namespace utopia {

    template <class FunctionSpace>
    class VelocityNewmarkIntegrator final : public NewmarkIntegrator<FunctionSpace>,
                                            public VelocityVariant<FunctionSpace> {
    public:
        using TimeDependentFunction = utopia::TimeDependentFunction<FunctionSpace>;
        using Super = utopia::NewmarkIntegrator<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

        template <class... Args>
        VelocityNewmarkIntegrator(Args &&...args) : Super(std::forward<Args>(args)...) {}

        const Vector_t &solution() const override { return this->x_old(); }

        void read(Input &in) override { Super::read(in); }

        bool update_IVP(const Vector_t &velocity) override {
            Vector_t x;
            update_x(velocity, x);
            return Super::update_IVP(x);
        }

        bool update_BVP() override {
            this->space()->apply_constraints(this->x_old());
            return true;
        }

        // The residual is the Newmark residual evaluated at x(v) = x_n + dt/2 (v_n + v).
        // The inertia term M a_{n+1} is added by NewmarkIntegrator::integrate_gradient, since
        // a_{n+1}(x(v)) = 4/dt^2 (x(v) - x_n - dt v_n) - a_n = 2/dt (v - v_n) - a_n
        bool hessian_and_gradient(const Vector_t &velocity, Matrix_t &H, Vector_t &g) const override {
            Vector_t x;
            update_x(velocity, x);
            return Super::hessian_and_gradient(x, H, g);
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

        bool value(const Vector_t &velocity, Scalar_t &v) const override {
            Vector_t x;
            update_x(velocity, x);
            return Super::value(x, v);
        }

        bool update(const Vector_t &velocity) override {
            Vector_t x;
            update_x(velocity, x);
            return Super::update(x);
        }

        ////////////////////////////////////////////////////////////////////////////////

        bool setup_IVP(Vector_t &x) override { return Super::setup_IVP(x); }

        void integrate_hessian(const Vector_t &x, Matrix_t &H) const override {
            UTOPIA_TRACE_REGION_BEGIN("VelocityNewmarkIntegrator::integrate_hessian");
            Super::integrate_hessian(x, H);
            H *= (this->delta_time() / 2);
            UTOPIA_TRACE_REGION_END("VelocityNewmarkIntegrator::integrate_hessian");
        }

        void hessian_from_velocity_to_displacement(Matrix_t &H) const { H *= (2 / this->delta_time()); }
        void hessian_from_displacement_to_velocity(Matrix_t &H) const { H *= (this->delta_time() / 2); }

        bool time_derivative(const Vector_t &x, Vector_t &dfdt) const override {
            return Super::time_derivative(x, dfdt);
        }

        bool report_solution(const Vector_t &) override { return Super::report_solution(solution()); }

        void initial_guess_for_solver(Vector_t &velocity) override {
            velocity.set(0.);
            this->impose_velocity_constraints(velocity);
        }

        void displacement(const Vector_t &velocity, Vector_t &result) const { update_x(velocity, result); }

    private:
        void update_x(const Vector_t &velocity, Vector_t &x) const { this->displacement_from_velocity(velocity, x); }
    };

}  // namespace utopia

#endif  // UTOPIA_VELOCITY_NEWMARK_INTEGRATOR_HPP
