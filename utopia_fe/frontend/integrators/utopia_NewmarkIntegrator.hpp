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

        bool setup_IVP(Vector_t &) override {
            if (!this->assemble_mass_matrix()) {
                return false;
            }

            this->space()->create_vector(active_stress_);
            active_stress_.set(0.);

            auto vlo = layout(active_stress_);

            x_old_.zeros(vlo);
            x_older_.zeros(vlo);
            internal_stress_old_.zeros(vlo);
            external_force_.zeros(vlo);
            return true;
        }

        bool update_IVP(const Vector_t &x) override {
            Super::update_IVP(x);

            const Scalar_t dt = this->delta_time();

            x_older_ = x_old_;
            x_old_ = x;

            active_stress_ = -internal_stress_old_ + (4. * external_force_);

            // In case of contact problems here we also have the contact stress
            this->gradient(x, internal_stress_old_);

            active_stress_ +=
                (4. / (dt * dt)) * ((*this->mass_matrix()) * (2. * x_old_ - x_older_)) - 2. * internal_stress_old_;

            return true;
        }

        template <class... Args>
        NewmarkIntegrator(Args &&... args) : Super(std::forward<Args>(args)...) {}

        virtual ~NewmarkIntegrator() = default;

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

        inline Vector_t &x_old() { return x_old_; }
        inline const Vector_t &x_old() const { return x_old_; }

        const Vector_t &solution() const override { return x_old(); }

    private:
        Vector_t x_old_, x_older_;
        Vector_t active_stress_, internal_stress_old_, external_force_;
    };

}  // namespace utopia

#endif  // UTOPIA_NEWMARK_INTEGRATOR_HPP
