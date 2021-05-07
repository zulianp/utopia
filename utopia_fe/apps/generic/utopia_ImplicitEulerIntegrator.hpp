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

        bool value(const Vector_t &x, Scalar_t &value) const override { return false; }

        bool gradient(const Vector_t &x, Vector_t &g) const override {
            if (!Super::gradient(x, g)) {
                return false;
            }

            integrate_gradient(x, g);
            return true;
        }

        bool update(const Vector_t &x) override { return true; }

        bool hessian(const Vector_t &x, Matrix_t &H) const override {
            if (!Super::hessian(x, H)) {
                return false;
            }

            integrate_hessian(x, H);
            return true;
        }

        bool hessian_and_gradient(const Vector_t &x, Matrix_t &H, Vector_t &g) const override {
            if (!Super::hessian_and_gradient(x, H, g)) {
                return false;
            }

            integrate_hessian(x, H);
            integrate_gradient(x, g);
            return true;
        }

        void read(Input &in) override {
            Super::read(in);

            // if (this->verbose()) {
            //     utopia::out() << "--------------------------------\n";
            //     this->describe(utopia::out().stream());
            //     utopia::out() << "--------------------------------\n";
            // }
        }

        bool setup_IVP(Vector_t &x) override { return this->assemble_mass_matrix(); }

        bool update_IVP(const Vector_t &) override { return true; }

        ImplicitEulerIntegrator(const std::shared_ptr<FunctionSpace> &space) : Super(space) {}

    private:
        void integrate_gradient(const Vector_t &x, Vector_t &g) const {
            const Scalar_t dt = this->delta_time();
            g *= dt;
            this->space()->apply_zero_constraints(g);
        }

        void integrate_hessian(const Vector_t &, Matrix_t &H) const {
            const Scalar_t dt = this->delta_time();

            system_ *= dt;
            system_ += (*this->mass_matrix());
        }
    };

}  // namespace utopia

#endif  // UTOPIA_IMPLICIT_EULER_INTEGRATOR_HPP
