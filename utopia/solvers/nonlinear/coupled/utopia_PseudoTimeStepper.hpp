#ifndef UTOPIA_PSEUDO_TIME_STEPPER_HPP
#define UTOPIA_PSEUDO_TIME_STEPPER_HPP

#include "utopia_Function.hpp"
#include "utopia_NewtonBase.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class TimeDependentFunction : public Function<Matrix, Vector> {
    public:
        using Traits = utopia::Traits<Matrix>;
        using Scalar = typename Traits::Scalar;

        using Function<Matrix, Vector>::update;

        virtual ~TimeDependentFunction() {}
        virtual bool update(const Scalar t) = 0;
    };

    template <class Matrix, class Vector>
    class PrototypeFunction : public TimeDependentFunction<Matrix, Vector> {
    public:
        using Traits = utopia::Traits<Matrix>;
        using Scalar = typename Traits::Scalar;

        PrototypeFunction(std::function<bool(const Vector &, Scalar &)> value,
                          std::function<bool(const Vector &, Vector &)> gradient,
                          std::function<bool(const Vector &, Matrix &)> hessian,
                          std::function<bool(const Scalar)> update,
                          std::function<void(Vector &x)> create_vector)
            : value_(value), gradient_(gradient), hessian_(hessian), update_(update), create_vector_(create_vector) {}

        bool hessian(const Vector &x, Matrix &H) const override { return hessian_(x, H); }
        bool value(const Vector &x, Scalar &value) const override { return value_(x, value); }
        bool gradient(const Vector &x, Vector &g) const override { return gradient_(x, g); }
        void create_vector(Vector &x) const { create_vector_(x); }
        bool update(const Scalar t) override { return update_(t); }

        std::function<bool(const Vector &, Scalar &)> value_;
        std::function<bool(const Vector &, Vector &)> gradient_;
        std::function<bool(const Vector &, Matrix &)> hessian_;
        std::function<bool(const Scalar)> update_;
        std::function<void(Vector &x)> create_vector_;
    };

    template <class Matrix, class Vector = typename Traits<Matrix>::Vector>
    class PseudoTimeStepper : public Configurable {
    public:
        using Traits = utopia::Traits<Matrix>;
        using SizeType = typename Traits::SizeType;
        using Scalar = typename Traits::Scalar;

        void read(Input &in) override {
            in.get("delta_time", delta_time_);
            in.get("start_time", start_time_);
            in.get("end_time", end_time_);
        }

        PseudoTimeStepper(const std::shared_ptr<NewtonInterface<Matrix, Vector>> &nonlinear_solver)
            : nonlinear_solver_(nonlinear_solver) {}

        void set_delta_time(const Scalar dt) { delta_time_ = dt; }
        void set_start_time(const Scalar t) { start_time_ = t; }
        void set_end_time(const Scalar t) { end_time_ = t; }

        bool solve(TimeDependentFunction<Matrix, Vector> &fun, Vector &x) {
            bool ok = false;
            for (Scalar t = start_time_; t <= end_time_; t += delta_time_) {
                fun.update(t);
                ok = nonlinear_solver_->solve(fun, x);
            }

            return ok;
        }

    private:
        std::shared_ptr<NewtonInterface<Matrix, Vector>> nonlinear_solver_;
        Scalar delta_time_{1};
        Scalar start_time_{0};
        Scalar end_time_{1};
    };
}  // namespace utopia

#endif  // UTOPIA_PSEUDO_TIME_STEPPER_HPP
