#ifndef UTOPIA_OBSTACLE_STABILIZED_NEWMARK_INTEGRATOR_HPP
#define UTOPIA_OBSTACLE_STABILIZED_NEWMARK_INTEGRATOR_HPP

#include "utopia_ContactFactory.hpp"
#include "utopia_ContactInterface.hpp"
#include "utopia_FEModelFunction.hpp"

#include <utility>

namespace utopia {

    template <class FunctionSpace>
    class ContactDependentFunction {
    public:
        virtual ~ContactDependentFunction() = default;
        virtual void set_contact(const std::shared_ptr<ContactInterface<FunctionSpace>> obstacle) = 0;
    };

    // Contact-stabilized Newmark method (gamma = 0.5, beta = 0.25), see
    // Deuflhard, Krause, Ertel, "A contact-stabilized Newmark method for dynamical contact problems", IJNME 2008.
    // The previous acceleration is replaced by the previous force, M a_n = -F(x_n), and the
    // predictor x_n + dt v_n is restricted to the admissible set.
    // Unlike Algorithm (21) of the paper, the velocity update uses the projected velocity
    // v_hat_n = (x_pred - x_n) / dt instead of v_n: v_{n+1} = v_hat_n + 2/dt (x_{n+1} - x_pred).
    // With v_n, the normal velocity of a node in persistent contact never changes (it keeps the impact
    // velocity while the node is at rest). With v_hat_n it vanishes, and the energy still does not increase,
    // since |v_hat_n|_M <= |v_n|_M (M-projection onto a convex set that contains 0)
    template <class FunctionSpace>
    class ObstacleStabilizedNewmark : public TimeDependentFunction<FunctionSpace>,
                                      public ContactDependentFunction<FunctionSpace> {
    public:
        using Super = utopia::TimeDependentFunction<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

        void read(Input &in) override {
            Super::read(in);

            in.get("debug", debug_);
            in.get("stabilized_formulation", stabilized_formulation_);
            in.get("consistent_initial_acceleration", consistent_initial_acceleration_);

            // if (stabilized_formulation_) {
            //     utopia::out() << "Using stabilized formulation!\n";
            // }
        }

        bool update_constraints(const Vector_t &x) {
            // utopia::out() << "ObstacleStabilizedNewmark::update_constraints\n";

            this->space()->displace(x);
            bool ok = obstacle_->assemble(*this->space());
            this->space()->displace(-x);
            return ok;
        }

        // Computes the predictor x_pred = x_n + dt v_hat_n, with the projected velocity v_hat_n
        void update_predictor() {
            const Scalar_t dt = this->delta_time();

            if (stabilized_formulation_) {
                assert(obstacle_);

                Vector_t temp;
                obstacle_->transform(velocity_old_, temp);
                temp = utopia::min(temp, (1 / dt) * obstacle_->gap());

                obstacle_->inverse_transform(temp, projected_velocity_);
            } else {
                projected_velocity_ = velocity_old_;
            }

            predictor_ = x_old_;
            predictor_ += dt * projected_velocity_;
        }

        bool setup_IVP(Vector_t &x) override {
            if (!this->assemble_mass_matrix()) {
                return false;
            }

            assert(this->mass_matrix());
            Scalar_t sum_mm = sum(*this->mass_matrix());
            has_zero_density_ = sum_mm == 0.0;

            auto vlo = layout(x);

            x_old_ = x;
            velocity_old_.zeros(vlo);
            force_old_.zeros(vlo);

            if (consistent_initial_acceleration_ && !has_zero_density_ && !compute_initial_force()) {
                return false;
            }

            if (stabilized_formulation_) {
                if (!obstacle_) {
                    Utopia::Abort(
                        "ObstacleStabilizedNewmark: stabilized_formulation requires a contact/obstacle (set_contact)!");
                }

                update_constraints(x);
            }

            update_predictor();
            return true;
        }

        // The force F(x_0) = r(x_0) - f(t_0) plays the role of -M a_0 (equation of motion at t_0)
        bool compute_initial_force() {
            if (this->has_time()) {
                this->space()->update(*this->time());
            }

            return this->function()->gradient(x_old_, force_old_);
        }

        // Incremental potential: E(x) = Psi(x) - f^T x + F_n^T x + 2/dt^2 (x - x_pred)^T M (x - x_pred)
        bool value(const Vector_t &x, Scalar_t &v) const override {
            if (!this->function()->value(x, v)) {
                return false;
            }

            if (!has_zero_density_) {
                const Scalar_t dt = this->delta_time();

                Vector_t dx = x - predictor_;
                Vector_t Mdx = (*this->mass_matrix()) * dx;
                const Scalar_t inertia = (2 / (dt * dt)) * dot(dx, Mdx);
                const Scalar_t previous_force = dot(force_old_, x);
                v += inertia + previous_force;
            }

            return true;
        }

        bool update_IVP(const Vector_t &x) override {
            Super::update_IVP(x);

            time_derivative(x, velocity_old_);

            this->function()->gradient(x, force_old_);

            // Store current solution
            x_old_ = x;

            if (stabilized_formulation_) {
                update_constraints(x);
            }

            update_predictor();
            return true;
        }

        bool time_derivative(const Vector_t &x, Vector_t &velocity) const override {
            const Scalar_t dt = this->delta_time();
            velocity = projected_velocity_;
            velocity += (2 / dt) * (x - predictor_);
            return true;
        }

        template <class... Args>
        ObstacleStabilizedNewmark(Args &&...args) : Super(std::forward<Args>(args)...) {}

        virtual ~ObstacleStabilizedNewmark() = default;

        void integrate_gradient(const Vector_t &x, Vector_t &g) const override {
            const Scalar_t dt2 = this->delta_time() * this->delta_time();

            if (!has_zero_density_) {
                Vector_t mom = (x - predictor_);
                mom *= (4.0 / dt2);
                g += (*this->mass_matrix()) * mom;
                g += force_old_;
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

        // The displacement is written by report_solution
        bool register_output(IO<FunctionSpace> &output) override {
            set_velocity_field();
            output.register_output_field(velocity_field_);
            return true;
        }

        bool update_output(IO<FunctionSpace> &output) override {
            set_velocity_field();
            output.update_output_field(velocity_field_);
            return true;
        }

        bool set_initial_condition(const Vector_t &x) override {
            x_old_ = x;
            return true;
        }

        inline void set_contact(const std::shared_ptr<ContactInterface<FunctionSpace>> obstacle) override {
            obstacle_ = obstacle;
        }

    protected:
        const Vector_t &velocity_old() const { return velocity_old_; }

    private:
        Vector_t x_old_, velocity_old_;
        Vector_t projected_velocity_;
        Vector_t predictor_;
        Field<FunctionSpace> velocity_field_{"velocity"};
        Vector_t force_old_;
        bool has_zero_density_{false};

        std::shared_ptr<ContactInterface<FunctionSpace>> obstacle_;
        bool debug_{false};
        bool stabilized_formulation_{true};
        bool consistent_initial_acceleration_{true};

        void set_velocity_field() {
            velocity_field_.set_space(this->space());
            velocity_field_.set_tensor_size(this->space()->n_var());
            velocity_field_.set_data(std::make_shared<Vector_t>(velocity_old_));
        }
    };

}  // namespace utopia

#endif  // UTOPIA_OBSTACLE_STABILIZED_NEWMARK_INTEGRATOR_HPP
