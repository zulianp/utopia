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

    // https://en.wikipedia.org/wiki/Newmark-beta_method
    // Unconditionally Stable gamma = 0.5, beta = 0.25
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

        void update_predictor() {
            assert(obstacle_);

            if (stabilized_formulation_) {
                const Scalar_t dt = this->delta_time();
                Vector_t temp;
                obstacle_->transform(velocity_old_, temp);
                temp *= dt;
                temp = utopia::min(temp, obstacle_->gap());

                obstacle_->inverse_transform(temp, predictor_);
            } else {
                predictor_ = this->delta_time() * velocity_old_;
            }

            predictor_ += x_old_;
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

            if (stabilized_formulation_) {
                update_constraints(x);
            }

            update_predictor();
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
            velocity = velocity_old_;
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
        Vector_t predictor_;
        Vector_t force_old_;
        bool has_zero_density_{false};

        std::shared_ptr<ContactInterface<FunctionSpace>> obstacle_;
        bool debug_{false};
        bool stabilized_formulation_{true};
    };

}  // namespace utopia

#endif  // UTOPIA_OBSTACLE_STABILIZED_NEWMARK_INTEGRATOR_HPP
