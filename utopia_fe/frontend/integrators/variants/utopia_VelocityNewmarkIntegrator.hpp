#ifndef UTOPIA_VELOCITY_NEWMARK_INTEGRATOR_HPP
#define UTOPIA_VELOCITY_NEWMARK_INTEGRATOR_HPP

#include "utopia_FEModelFunction.hpp"
#include "utopia_NewmarkIntegrator.hpp"
#include "utopia_VelocityVariant.hpp"

#include <utility>

namespace utopia {

    template <class FunctionSpace>
    class VelocityNewmarkIntegrator final : public NewmarkIntegrator<FunctionSpace>,
                                            public VelocityVariant<FunctionSpace> {
    public:
        using Super = utopia::NewmarkIntegrator<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

        template <class... Args>
        VelocityNewmarkIntegrator(Args &&... args) : Super(std::forward<Args>(args)...) {}

        const Vector_t &solution() const override { return this->x_old(); }

        void read(Input &in) override { Super::read(in); }

        bool update_IVP(const Vector_t &velocity) override {
            Vector_t x = this->x_old();
            update_x(velocity, x);
            velocity_old_ = velocity;
            return Super::update_IVP(x);
        }

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

        ////////////////////////////////////////////////////////////////////////////////

        bool setup_IVP(Vector_t &x) override {
            velocity_old_.zeros(layout(x));
            return Super::setup_IVP(x);
        }

        void integrate_gradient(const Vector_t &x, Vector_t &g) const override {
            Super::integrate_gradient(x, g);
            g *= 2. / this->delta_time();
        }

        void integrate_hessian(const Vector_t &x, Matrix_t &H) const override { Super::integrate_hessian(x, H); }

        bool time_derivative(const Vector_t &x, Vector_t &dfdt) const override {
            return Super::time_derivative(x, dfdt);
        }

        bool report_solution(const Vector_t &) override { return Super::report_solution(solution()); }

    private:
        void update_x(const Vector_t &velocity, Vector_t &x) const {
            x = this->x_old();
            x += (0.5 * this->delta_time()) * (velocity_old_ + velocity);
        }

        Vector_t velocity_old_;
    };

}  // namespace utopia

#endif  // UTOPIA_VELOCITY_NEWMARK_INTEGRATOR_HPP
