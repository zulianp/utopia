#ifndef UTOPIA_NEWMARK_INTEGRATOR_HPP
#define UTOPIA_NEWMARK_INTEGRATOR_HPP

#include "utopia_FEModelFunction.hpp"
#include "utopia_polymorphic_LinearSolver.hpp"

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
        using Layout_t = typename Traits<FunctionSpace>::Layout;

        class State {
        public:
            Field<FunctionSpace> x{"disp"}, velocity{"velocity"}, acceleration{"acceleration"};
            bool has_zero_density{false};
        };

        void read(Input &in) override {
            Super::read(in);
            in.get("enable_restart", enable_restart_);
            if (enable_restart_) {
                export_time_integrator_fields_ = enable_restart_;
            }

            if (!export_time_integrator_fields_) {
                in.get("output",
                       [&](Input &node) { node.get("export_time_integrator_fields", export_time_integrator_fields_); });
            }

            in.get("fail_on_field_not_found", fail_on_field_not_found_);
            in.get("consistent_initial_acceleration", consistent_initial_acceleration_);
        }

        bool setup_IVP(IO<FunctionSpace> &input) override {
            if (!this->assemble_mass_matrix()) {
                return false;
            }

            assert(this->mass_matrix());
            Scalar_t sum_mm = sum(*this->mass_matrix());
            state_->has_zero_density = sum_mm == 0.0;

            state_->x.set_space(this->space());
            state_->velocity.set_space(this->space());
            state_->acceleration.set_space(this->space());

            bool ok = input.read_nodal(state_->x, fail_on_field_not_found_) &&
                      input.read_nodal(state_->velocity, fail_on_field_not_found_) &&
                      input.read_nodal(state_->acceleration, fail_on_field_not_found_);

            if (!ok) {
                auto vlo = layout(this->space()->comm(), this->space()->n_local_dofs(), this->space()->n_dofs());
                init_vectors(vlo);
            }

            return ok;
        }

        void init_vectors(const Layout_t &vlo) {
            if (state_->x.empty()) {
                state_->x.set_space(this->space());
                state_->x.set_tensor_size(this->space()->n_var());
                state_->x.set_data(std::make_shared<Vector_t>(vlo, 0));
            }

            if (state_->velocity.empty()) {
                state_->velocity.set_space(this->space());
                state_->velocity.set_tensor_size(this->space()->n_var());
                state_->velocity.set_data(std::make_shared<Vector_t>(vlo, 0));
            }

            if (state_->acceleration.empty()) {
                state_->acceleration.set_space(this->space());
                state_->acceleration.set_tensor_size(this->space()->n_var());
                state_->acceleration.set_data(std::make_shared<Vector_t>(vlo, 0));
            }
        }

        bool register_output(IO<FunctionSpace> &output) override {
            output.register_output_field(state_->x);

            // We do not need all of this on disk if no restart is required
            if (export_time_integrator_fields_) {
                output.register_output_field(state_->velocity);
                output.register_output_field(state_->acceleration);
            }

            return true;
        }

        bool update_output(IO<FunctionSpace> &output) override {
            output.update_output_field(state_->x);

            // We do not need all of this on disk if no restart is required
            if (export_time_integrator_fields_) {
                output.update_output_field(state_->velocity);
                output.update_output_field(state_->acceleration);
            }

            return true;
        }

        bool setup_IVP(Vector_t &x) override {
            if (!this->assemble_mass_matrix()) {
                return false;
            }

            assert(this->mass_matrix());
            Scalar_t sum_mm = sum(*this->mass_matrix());
            state_->has_zero_density = sum_mm == 0.0;

            auto vlo = layout(x);
            // init_vectors(vlo);

            state_->x.set_space(this->space());
            state_->velocity.set_space(this->space());
            state_->acceleration.set_space(this->space());

            state_->x.set_tensor_size(this->space()->n_var());
            state_->velocity.set_tensor_size(this->space()->n_var());
            state_->acceleration.set_tensor_size(this->space()->n_var());

            state_->x.set_data(std::make_shared<Vector_t>(x));
            state_->velocity.set_data(std::make_shared<Vector_t>(vlo, 0));
            state_->acceleration.set_data(std::make_shared<Vector_t>(vlo, 0));

            if (consistent_initial_acceleration_ && !state_->has_zero_density) {
                return compute_initial_acceleration();
            }

            return true;
        }

        // Solves M a_0 = f(t_0) - r(x_0), i.e., the equation of motion at t_0.
        // On the Dirichlet dofs a_0 = 0.
        bool compute_initial_acceleration() {
            if (this->has_time()) {
                this->space()->update(*this->time());
            }

            Vector_t rhs;
            if (!this->function()->gradient(x_old(), rhs)) {
                return false;
            }

            rhs *= -1.0;
            this->space()->apply_zero_constraints(rhs);

            Matrix_t M = *this->mass_matrix();
            this->space()->apply_constraints(M);

            OmniLinearSolver<Matrix_t, Vector_t> solver;
            acceleration_old().set(0.);
            return solver.solve(M, rhs, acceleration_old());
        }

        // Incremental potential of the Newmark step:
        // E(x) = Psi(x) - f^T x + 2/dt^2 (x - x_pred)^T M (x - x_pred), x_pred = x_n + dt v_n + dt^2/4 a_n
        bool value(const Vector_t &x, Scalar_t &v) const override {
            if (!this->function()->value(x, v)) {
                return false;
            }

            if (!state_->has_zero_density) {
                const Scalar_t dt = this->delta_time();
                const Scalar_t dt2 = dt * dt;

                Vector_t dx = x - x_old();
                dx -= dt * velocity_old();
                dx -= (dt2 / 4) * acceleration_old();

                Vector_t Mdx = (*this->mass_matrix()) * dx;
                v += (2 / dt2) * dot(dx, Mdx);
            }

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
            this->space()->apply_constraints_time_derivative(this->velocity_old());
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

        inline Vector_t &x_old() { return state_->x.data(); }
        inline const Vector_t &x_old() const { return state_->x.data(); }

        const Vector_t &solution() const override { return x_old(); }
        const Vector_t &velocity() const { return state_->velocity.data(); }
        const Vector_t &acceleration() const { return state_->acceleration.data(); }

        Vector_t &velocity() { return state_->velocity.data(); }
        Vector_t &acceleration() { return state_->acceleration.data(); }

        bool set_initial_condition(const Vector_t &x) override {
            x_old() = x;
            return true;
        }

        inline std::shared_ptr<State> state() { return state_; }
        inline void set_state(const std::shared_ptr<State> &state) { state_ = state; }

    protected:
        const Vector_t &velocity_old() const { return velocity(); }
        Vector_t &velocity_old() { return velocity(); }
        const Vector_t &acceleration_old() const { return acceleration(); }
        Vector_t &acceleration_old() { return acceleration(); }

        // Trapezoidal rule x_{n+1} = x_n + dt/2 (v_n + v_{n+1}), used by the velocity variants
        void displacement_from_velocity(const Vector_t &velocity, Vector_t &x) const {
            x = x_old();
            x += (0.5 * this->delta_time()) * (velocity_old() + velocity);
        }

        // Sets the Dirichlet dofs of the velocity to v_{n+1} = 2/dt (g(t_{n+1}) - x_n) - v_n,
        // so that the displacement obtained with displacement_from_velocity satisfies x_{n+1} = g(t_{n+1})
        void impose_velocity_constraints(Vector_t &velocity) const {
            Vector_t g = x_old();
            this->space()->apply_constraints(g);

            Vector_t v_bc = (2 / this->delta_time()) * (g - x_old());
            v_bc -= velocity_old();

            // 1 on the free dofs and 0 on the Dirichlet dofs
            Vector_t is_free(layout(velocity), 1.);
            this->space()->apply_zero_constraints(is_free);

            Vector_t v_free = e_mul(is_free, velocity);
            Vector_t v_bc_free = e_mul(is_free, v_bc);

            velocity = v_free;
            velocity += v_bc;
            velocity -= v_bc_free;
        }

    private:
        bool enable_restart_{false};
        bool export_time_integrator_fields_{true};
        bool fail_on_field_not_found_{true};
        bool consistent_initial_acceleration_{true};
        std::shared_ptr<State> state_;
    };

}  // namespace utopia

#endif  // UTOPIA_NEWMARK_INTEGRATOR_HPP
