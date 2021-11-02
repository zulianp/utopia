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

            assert(this->mass_matrix());
            Scalar_t sum_mm = sum(*this->mass_matrix());
            has_zero_density_ = sum_mm == 0.0;

            auto vlo = layout(x);

            // x_old_.zeros(vlo);
            x_old_ = x;
            velocity_old_.zeros(vlo);
            acceleration_old_.zeros(vlo);
            return true;
        }

        bool update_IVP(const Vector_t &x) override {
            Super::update_IVP(x);

            time_second_derivative(x, acceleration_old_);
            time_derivative(x, velocity_old_);

            // Store current solution
            x_old_ = x;
            return true;
        }

        bool time_second_derivative(const Vector_t &x, Vector_t &acceleration) const {
            const Scalar_t dt = this->delta_time();
            const Scalar_t dt2 = dt * dt;

            acceleration = -acceleration_old_;
            acceleration += (4 / dt2) * (x - x_old_ - dt * velocity_old_);
            return true;
        }

        bool time_derivative(const Vector_t &x, Vector_t &velocity) const override {
            velocity = -velocity_old_;
            velocity += (2 / this->delta_time()) * (x - x_old_);
            return true;
        }

        template <class... Args>
        NewmarkIntegrator(Args &&... args) : Super(std::forward<Args>(args)...) {}

        virtual ~NewmarkIntegrator() = default;

        void integrate_gradient(const Vector_t &x, Vector_t &g) const override {
            const Scalar_t dt2 = this->delta_time() * this->delta_time();

            if (!has_zero_density_) {
                Vector_t mom = (x - x_old_);
                mom -= this->delta_time() * velocity_old_;
                mom *= (4.0 / dt2);
                mom -= acceleration_old_;

                g += (*this->mass_matrix()) * mom;
            }

            this->space()->apply_zero_constraints(g);
        }

        void integrate_hessian(const Vector_t &, Matrix_t &H) const override {
            if (!has_zero_density_) {
                const Scalar_t dt2 = this->delta_time() * this->delta_time();
                H += (4. / dt2) * (*this->mass_matrix());
            }

            this->space()->apply_constraints(H);
        }

        inline Vector_t &x_old() { return x_old_; }
        inline const Vector_t &x_old() const { return x_old_; }

        const Vector_t &solution() const override { return x_old(); }
        const Vector_t &velocity() const { return velocity_old_; }
        const Vector_t &acceleration() const { return acceleration_old_; }

        bool set_initial_condition(const Vector_t &x) override {
            x_old_ = x;
            return true;
        }

    protected:
        const Vector_t &velocity_old() const { return velocity_old_; }

    private:
        Vector_t x_old_, velocity_old_, acceleration_old_;
        bool has_zero_density_{false};
    };

}  // namespace utopia

#endif  // UTOPIA_NEWMARK_INTEGRATOR_HPP
