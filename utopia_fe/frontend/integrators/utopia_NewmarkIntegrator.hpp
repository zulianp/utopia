#ifndef UTOPIA_NEWMARK_INTEGRATOR_HPP
#define UTOPIA_NEWMARK_INTEGRATOR_HPP

#include "utopia_FEModelFunction.hpp"

#include <utility>

namespace utopia {

    // https://en.wikipedia.org/wiki/Newmark-beta_method
    // Unconditionally Stable gamma = 0.5, beta = 0.25
    template <class FunctionSpace>
    class NewmarkIntegrator : public TimeDependentFunction<FunctionSpace> {
    public:
        using Super = utopia::TimeDependentFunction<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

        void read(Input &in) override { Super::read(in); }

        bool setup_IVP(Vector_t &x) override {
            if (!this->assemble_mass_matrix()) {
                return false;
            }

            auto vlo = layout(x);

            x_old_.zeros(vlo);
            velocity_old_.zeros(vlo);
            acceleration_old_.zeros(vlo);
            return true;
        }

        bool update_IVP(const Vector_t &x) override {
            Super::update_IVP(x);

            const Scalar_t dt = this->delta_time();

            Vector_t acceleration_new;
            this->gradient(x, acceleration_new);
            acceleration_new *= -1.0;

            velocity_old_ += dt / 2.0 * (acceleration_new + acceleration_old_);
            x_old_ = x;
            acceleration_old_ = acceleration_new;
            return true;
        }

        template <class... Args>
        NewmarkIntegrator(Args &&... args) : Super(std::forward<Args>(args)...) {}

        virtual ~NewmarkIntegrator() = default;

        void integrate_gradient(const Vector_t &x, Vector_t &g) const override {
            const Scalar_t dt2 = this->delta_time() * this->delta_time();

            Vector_t mom = (*this->mass_matrix()) * (x - x_old_ - this->delta_time() * velocity_old_);
            g -= acceleration_old_;
            g *= dt2 / 4.0;
            g += mom;

            this->space()->apply_zero_constraints(g);
        }

        void integrate_hessian(const Vector_t &, Matrix_t &H) const override {
            const Scalar_t dt2 = this->delta_time() * this->delta_time();
            H *= (dt2 / 4.);
            H += (*this->mass_matrix());
            this->space()->apply_constraints(H);
        }

        inline Vector_t &x_old() { return x_old_; }
        inline const Vector_t &x_old() const { return x_old_; }

        const Vector_t &solution() const override { return x_old(); }

    private:
        Vector_t x_old_, velocity_old_, acceleration_old_;
    };

}  // namespace utopia

#endif  // UTOPIA_NEWMARK_INTEGRATOR_HPP
