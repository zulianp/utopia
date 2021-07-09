#ifndef UTOPIA_IMPLICIT_EULER_INTEGRATOR_HPP
#define UTOPIA_IMPLICIT_EULER_INTEGRATOR_HPP

#include "utopia_FEModelFunction.hpp"

#include <utility>

namespace utopia {

    template <class FunctionSpace>
    class ImplicitEulerIntegrator : public TimeDependentFunction<FunctionSpace> {
    public:
        using Super = utopia::TimeDependentFunction<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

        template <class... Args>
        ImplicitEulerIntegrator(Args &&... args) : Super(std::forward<Args>(args)...) {}

        void read(Input &in) override { Super::read(in); }

        bool setup_IVP(Vector_t &x) override {
            if (!this->assemble_mass_matrix()) {
                return false;
            }
            x_old_.zeros(layout(x));
            assert(!empty(x_old_));
            return true;
        }

        bool update_IVP(const Vector_t &x) override {
            Super::update_IVP(x);
            x_old_ = x;
            assert(!empty(x_old_));
            return true;
        }

        void integrate_gradient(const Vector_t &x, Vector_t &g) const override {
            const Scalar_t dt = this->delta_time();

            assert(!empty(g));
            assert(!empty(x_old_));

            g *= dt;
            g += (*this->mass_matrix()) * (x - x_old_);
            this->space()->apply_zero_constraints(g);
        }

        bool time_derivative(const Vector_t &x, Vector_t &dfdt) const override {
            dfdt = (x - x_old_);
            dfdt *= 1. / this->delta_time();
            return true;
        }

        void integrate_hessian(const Vector_t &, Matrix_t &H) const override {
            const Scalar_t dt = this->delta_time();
            H *= dt;
            H += (*this->mass_matrix());
        }

        inline Vector_t &x_old() { return x_old_; }
        inline const Vector_t &x_old() const { return x_old_; }

    private:
        Vector_t x_old_;
    };

}  // namespace utopia

#endif  // UTOPIA_IMPLICIT_EULER_INTEGRATOR_HPP
