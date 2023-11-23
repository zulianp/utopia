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

        virtual ~TimeDependentFunction() {}
        virtual bool update(const Scalar t) = 0;
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

        PseudoTimeStepper(const std::shared_ptr<NewtonBase<Matrix, Vector>> &nonlinear_solver)
            : nonlinear_solver_(nonlinear_solver) {}

        void set_delta_time(const Scalar dt) { delta_time_ = dt; }
        void set_start_time(const Scalar t) { start_time_ = t; }
        void set_end_time(const Scalar t) { end_time_ = t; }

        bool solve(TimeDependentFunction<Matrix, Vector> &fun, Vector &x) {
            bool ok = false;
            for (Scalar t = start_time_; t < end_time_; t += delta_time_) {
                fun.update(t);
                ok = nonlinear_solver_->solve(fun, x);
            }

            return ok;
        }

    private:
        std::shared_ptr<NewtonBase<Matrix, Vector>> nonlinear_solver_;
        Scalar delta_time_{1};
        Scalar start_time_{0};
        Scalar end_time_{1};
    };
}  // namespace utopia

#endif  // UTOPIA_PSEUDO_TIME_STEPPER_HPP
