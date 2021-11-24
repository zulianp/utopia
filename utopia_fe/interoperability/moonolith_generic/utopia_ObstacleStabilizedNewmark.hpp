#ifndef UTOPIA_OBSTACLE_STABILIZED_NEWMARK_INTEGRATOR_HPP
#define UTOPIA_OBSTACLE_STABILIZED_NEWMARK_INTEGRATOR_HPP

#include "utopia_FEModelFunction.hpp"
#include "utopia_IObstacle.hpp"
#include "utopia_ObstacleFactory.hpp"

#include <utility>

namespace utopia {

    template <class FunctionSpace>
    class ObstacleDependentFunction {
    public:
        virtual ~ObstacleDependentFunction() = default;
        virtual void set_obstacle(const std::shared_ptr<IObstacle<FunctionSpace>> obstacle) = 0;
    };

    // https://en.wikipedia.org/wiki/Newmark-beta_method
    // Unconditionally Stable gamma = 0.5, beta = 0.25
    template <class FunctionSpace>
    class ObstacleStabilizedNewmark : public TimeDependentFunction<FunctionSpace>,
                                      public ObstacleDependentFunction<FunctionSpace> {
    public:
        using Super = utopia::TimeDependentFunction<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

        void read(Input &in) override {
            Super::read(in);

            // in.get("linear_obstacle", linear_obstacle_);
            in.get("debug", debug_);

            // Obstacle is initialized outside
            // {
            //     if (!obstacle_) {
            //         std::string type;
            //         in.get("obstacle",
            //                [&](Input &node) { obstacle_ = ObstacleFactory<FunctionSpace>::new_obstacle(node); });
            //     }
            // }
        }

        bool update_constraints(const Vector_t &x) {
            utopia::out() << "update_constraints\n";

            this->space()->displace(x);
            bool ok = obstacle_->assemble(*this->space());

            // if (debug_) {
            //     static int iter_debug = 0;
            //     ouput_debug_data(iter_debug++, x);
            // }

            this->space()->displace(-x);
            return ok;
        }

        void update_predictor() {
            assert(obstacle_);

            const Scalar_t dt = this->delta_time();
            Vector_t temp;
            obstacle_->transform(velocity_old_, temp);
            temp *= dt;
            temp = utopia::min(temp, obstacle_->gap());

            obstacle_->inverse_transform(temp, delta_predictor_);
        }

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

            update_constraints(x);
            update_predictor();
            return true;
        }

        bool update_IVP(const Vector_t &x) override {
            Super::update_IVP(x);

            time_second_derivative(x, acceleration_old_);
            time_derivative(x, velocity_old_);

            // Store current solution
            x_old_ = x;

            update_constraints(x);
            update_predictor();
            return true;
        }

        bool time_second_derivative(const Vector_t &x, Vector_t &acceleration) const {
            const Scalar_t dt = this->delta_time();
            const Scalar_t dt2 = dt * dt;

            acceleration = -acceleration_old_;
            acceleration += (4 / dt2) * (x - x_old_ - delta_predictor_);
            return true;
        }

        bool time_derivative(const Vector_t &x, Vector_t &velocity) const override {
            velocity = -velocity_old_;
            velocity += (2 / this->delta_time()) * (x - x_old_ - delta_predictor_);
            return true;
        }

        bool time_second_derivative_contact(const Vector_t &x, const Vector_t &acceleration, Vector_t &a_con) {
            const Scalar_t dt = this->delta_time();
            a_con = 2 / dt * (x_old_ + delta_predictor_ - x);
            a_con += 0.5 * (acceleration + acceleration_old_);
        }

        template <class... Args>
        ObstacleStabilizedNewmark(Args &&... args) : Super(std::forward<Args>(args)...) {}

        virtual ~ObstacleStabilizedNewmark() = default;

        void integrate_gradient(const Vector_t &x, Vector_t &g) const override {
            const Scalar_t dt2 = this->delta_time() * this->delta_time();

            if (!has_zero_density_) {
                Vector_t mom = (x - x_old_);
                mom -= delta_predictor_;
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

        inline void set_obstacle(const std::shared_ptr<IObstacle<FunctionSpace>> obstacle) override {
            obstacle_ = obstacle;
        }

    protected:
        const Vector_t &velocity_old() const { return velocity_old_; }

    private:
        Vector_t x_old_, velocity_old_, acceleration_old_;
        Vector_t delta_predictor_;
        bool has_zero_density_{false};

        std::shared_ptr<IObstacle<FunctionSpace>> obstacle_;
        bool debug_{false};
    };

}  // namespace utopia

#endif  // UTOPIA_OBSTACLE_STABILIZED_NEWMARK_INTEGRATOR_HPP
