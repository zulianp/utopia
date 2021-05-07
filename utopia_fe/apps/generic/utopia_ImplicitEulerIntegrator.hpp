#ifndef UTOPIA_IMPLICIT_EULER_INTEGRATOR_HPP
#define UTOPIA_IMPLICIT_EULER_INTEGRATOR_HPP

#include "utopia_FEModelFunction.hpp"

namespace utopia {

    template <class FunctionSpace>
    class ImplicitEulerIntegrator final : public TimeDependentFunction<FunctionSpace> {
    public:
        using Super = utopia::TimeDependentFunction<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

        ImplicitEulerIntegrator(const std::shared_ptr<FunctionSpace> &space) : Super(space) {}

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
            mass_times_x_old_.zeros(layout(x));
            assert(!empty(mass_times_x_old_));
            return true;
        }

        bool update_IVP(const Vector_t &x) override {
            mass_times_x_old_ = (*this->mass_matrix()) * x;
            assert(!empty(mass_times_x_old_));
            return true;
        }

        void integrate_gradient(const Vector_t &x, Vector_t &g) const override {
            const Scalar_t dt = this->delta_time();

            assert(!empty(g));
            assert(!empty(mass_times_x_old_));

            g *= dt;
            g -= mass_times_x_old_;
            g += (*this->mass_matrix()) * x;
            this->space()->apply_zero_constraints(g);
        }

        void integrate_hessian(const Vector_t &, Matrix_t &H) const override {
            const Scalar_t dt = this->delta_time();
            H *= dt;
            H += (*this->mass_matrix());
        }

    private:
        Vector_t mass_times_x_old_;
    };

}  // namespace utopia

#endif  // UTOPIA_IMPLICIT_EULER_INTEGRATOR_HPP
