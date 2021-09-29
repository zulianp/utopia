#ifndef UTOPIA_INCREMENTAL_LOADING_HPP
#define UTOPIA_INCREMENTAL_LOADING_HPP

#include "utopia_FEModelFunction.hpp"

#include <utility>

namespace utopia {

    template <class FunctionSpace>
    class IncrementalLoading : public TimeDependentFunction<FunctionSpace> {
    public:
        using Super = utopia::TimeDependentFunction<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

        template <class... Args>
        IncrementalLoading(Args &&... args) : Super(std::forward<Args>(args)...) {}

        virtual ~IncrementalLoading() = default;

        void read(Input &in) override { Super::read(in); }

        bool setup_IVP(Vector_t &x) override {
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

        void integrate_gradient(const Vector_t &, Vector_t &g) const override {
            this->space()->apply_zero_constraints(g);
        }

        bool time_derivative(const Vector_t &x, Vector_t &dfdt) const override {
            dfdt = (x - x_old_);
            dfdt *= 1. / this->delta_time();
            return true;
        }

        void integrate_hessian(const Vector_t &, Matrix_t &) const override {}

        inline Vector_t &x_old() { return x_old_; }
        inline const Vector_t &x_old() const { return x_old_; }
        const Vector_t &solution() const override { return x_old(); }

        bool set_initial_condition(const Vector_t &x) override {
            x_old_ = x;
            return true;
        }

    private:
        Vector_t x_old_;
    };

}  // namespace utopia

#endif  // UTOPIA_INCREMENTAL_LOADING_HPP
