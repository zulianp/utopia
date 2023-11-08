#ifndef UTOPIA_NEWMARK_INTEGRATOR_HPP
#define UTOPIA_NEWMARK_INTEGRATOR_HPP

#include "utopia_FEModelFunction.hpp"

#include <memory>
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

        class State {
        public:
            Vector_t x, velocity, acceleration;
            bool has_zero_density{false};
        };

        void read(Input &in) override { Super::read(in); }

        bool setup_IVP(Vector_t &x) override {
            if (!this->assemble_mass_matrix()) {
                return false;
            }

            assert(this->mass_matrix());
            Scalar_t sum_mm = sum(*this->mass_matrix());
            state_->has_zero_density = sum_mm == 0.0;

            auto vlo = layout(x);

            // x_old_.zeros(vlo);
            x_old() = x;
            velocity_old().zeros(vlo);
            acceleration_old().zeros(vlo);
            return true;
        }

        bool update_IVP(const Vector_t &x) override {
            Super::update_IVP(x);

            time_second_derivative(x, acceleration_old());
            time_derivative(x, velocity_old());

            // Store current solution
            x_old() = x;
            return true;
        }

        bool update_BVP() override {
            this->space()->apply_constraints(this->x_old());
            return true;
        }

        bool time_second_derivative(const Vector_t &x, Vector_t &acceleration) const {
            const Scalar_t dt = this->delta_time();
            const Scalar_t dt2 = dt * dt;

            acceleration = -acceleration_old();
            acceleration += (4 / dt2) * (x - x_old() - dt * velocity_old());
            return true;
        }

        bool time_derivative(const Vector_t &x, Vector_t &velocity) const override {
            velocity = -velocity_old();
            velocity += (2 / this->delta_time()) * (x - x_old());
            return true;
        }

        template <class... Args>
        NewmarkIntegrator(Args &&...args) : Super(std::forward<Args>(args)...), state_(std::make_shared<State>()) {}

        virtual ~NewmarkIntegrator() = default;

        void integrate_gradient(const Vector_t &x, Vector_t &g) const override {
            const Scalar_t dt2 = this->delta_time() * this->delta_time();

            if (!state_->has_zero_density) {
                Vector_t mom = (x - x_old());
                mom -= this->delta_time() * velocity_old();
                mom *= (4.0 / dt2);
                mom -= acceleration_old();

                g += (*this->mass_matrix()) * mom;
            }

            this->space()->apply_zero_constraints(g);
        }

        void integrate_hessian(const Vector_t &, Matrix_t &H) const override {
            if (!state_->has_zero_density) {
                const Scalar_t dt2 = this->delta_time() * this->delta_time();
                H += (4. / dt2) * (*this->mass_matrix());
            }

            this->space()->apply_constraints(H);
        }

        inline Vector_t &x_old() { return state_->x; }
        inline const Vector_t &x_old() const { return state_->x; }

        const Vector_t &solution() const override { return x_old(); }
        const Vector_t &velocity() const { return state_->velocity; }
        const Vector_t &acceleration() const { return state_->acceleration; }

        bool set_initial_condition(const Vector_t &x) override {
            x_old() = x;
            return true;
        }

        inline std::shared_ptr<State> state() { return state_; }
        inline void set_state(const std::shared_ptr<State> &state) { state_ = state; }

    protected:
        const Vector_t &velocity_old() const { return state_->velocity; }
        Vector_t &velocity_old() { return state_->velocity; }
        const Vector_t &acceleration_old() const { return state_->acceleration; }
        Vector_t &acceleration_old() { return state_->acceleration; }

    private:
        std::shared_ptr<State> state_;
    };

}  // namespace utopia

#endif  // UTOPIA_NEWMARK_INTEGRATOR_HPP
