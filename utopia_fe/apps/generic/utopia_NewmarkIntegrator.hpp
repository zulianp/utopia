#ifndef UTOPIA_NEWMARK_INTEGRATOR_HPP
#define UTOPIA_NEWMARK_INTEGRATOR_HPP

#include "utopia_FEModelFunction.hpp"

namespace utopia {

    // https://en.wikipedia.org/wiki/Newmark-beta_method
    // Unconditionally Stable gamma = 0.5, beta = 0.25
    template <class FunctionSpace>
    class NewmarkIntegrator final : public TimeDependentFunction<FunctionSpace> {
    public:
        using Super = utopia::TimeDependentFunction<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

        void read(Input &in) override {
            Super::read(in);

            // if (this->verbose()) {
            //     utopia::out() << "--------------------------------\n";
            //     this->describe(utopia::out().stream());
            //     utopia::out() << "--------------------------------\n";
            // }
        }

        bool setup_IVP(Vector_t &x) override {
            if (!this->assemble_mass_matrix()) {
                return false;
            }

            this->space()->create_vector(active_stress_);
            active_stress_.set(0.);

            auto vlo = layout(active_stress_);

            x_old_.zeros(vlo);
            x_older_.zeros(vlo);
            internal_stress_old_.zeros(vlo);
            internal_stress_older_.zeros(vlo);
            external_force_.zeros(vlo);
            internal_stress_.zeros(vlo);
            return true;
        }

        bool update_IVP(const Vector_t &x) override {
            const Scalar_t dt = this->delta_time();

            x_older_ = x_old_;
            x_old_ = x;

            // this->function()->gradient(x, internal_stress_);
            this->gradient(x, internal_stress_);

            internal_stress_older_ = internal_stress_old_;
            internal_stress_old_ = internal_stress_;

            active_stress_ = (4. / (dt * dt)) * ((*this->mass_matrix()) * (2. * x_old_ - x_older_)) -
                             2. * internal_stress_old_ - internal_stress_older_ + (4. * external_force_);

            return true;
        }

        NewmarkIntegrator(const std::shared_ptr<FunctionSpace> &space) : Super(space) {}

        void integrate_gradient(const Vector_t &x, Vector_t &g) const override {
            const Scalar_t dt2 = this->delta_time() * this->delta_time();

            g -= active_stress_;
            g *= (dt2 / 4.);
            g += ((*this->mass_matrix()) * x);
            this->space()->apply_zero_constraints(g);
        }

        void integrate_hessian(const Vector_t &, Matrix_t &H) const override {
            const Scalar_t dt2 = this->delta_time() * this->delta_time();
            H *= (dt2 / 4.);
            H += (*this->mass_matrix());
            this->space()->apply_constraints(H);
        }

    private:
        Vector_t x_old_, x_older_;
        Vector_t active_stress_, internal_stress_, internal_stress_old_, internal_stress_older_, external_force_;
        // Scalar_t beta_{0.25};
        // Scalar_t gamma_{0.5};
    };

}  // namespace utopia

#endif  // UTOPIA_NEWMARK_INTEGRATOR_HPP