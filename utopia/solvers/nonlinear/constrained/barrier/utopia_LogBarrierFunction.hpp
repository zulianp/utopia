#ifndef UTOPIA_LOG_BARRIER_FUNCTION_HPP
#define UTOPIA_LOG_BARRIER_FUNCTION_HPP

#include "utopia_BoxConstraints.hpp"
#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_Layout.hpp"
#include "utopia_Options.hpp"

#include <limits>

namespace utopia {
    template <class Matrix, class Vector>
    class LogBarrierFunction : public Function<Matrix, Vector> {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Function = utopia::Function<Matrix, Vector>;
        using BoxConstraints = utopia::BoxConstraints<Vector>;
        using Super = Function;

        void read(Input &in) override {
            if (!Options()
                     .add_option(
                         "barrier_parameter", barrier_parameter_, "see Numerical Optimization - J. Nocedal, S. Wright.")
                     .add_option("barrier_parameter_shrinking_factor",
                                 barrier_parameter_shrinking_factor_,
                                 "Factor with which the barrier term is reduced.")
                     .add_option("min_barrier_parameter", min_barrier_parameter_, "Smallest barrier parameter allowed.")
                     .add_option(
                         "soft_boundary", soft_boundary_, "A value in (0,0.01) used to project in the feasible region.")
                     .add_option("verbose", verbose_, "Enable/Disable verbose output.")
                     .parse(in)) {
                return;
            }

            reset();
        }

        LogBarrierFunction() {}

        inline void set_unconstrained_function(const std::shared_ptr<Function> &unconstrained) {
            unconstrained_ = unconstrained;
        }

        void set_box_constraints(const std::shared_ptr<BoxConstraints> &box) { box_ = box; }

        LogBarrierFunction(const std::shared_ptr<Function> &unconstrained, const std::shared_ptr<BoxConstraints> &box)
            : unconstrained_(unconstrained), box_(box) {}

        bool hessian(const Vector &x, Matrix &H) const override {
            if (!unconstrained_->hessian(x, H)) {
                return false;
            }

            extend_hessian(x, H);
            return true;
        }
        bool hessian(const Vector &x, Matrix &H, Matrix &preconditioner) const override {
            if (!unconstrained_->hessian(x, H, preconditioner)) {
                return false;
            }

            extend_hessian(x, H);
            return true;
        }

        bool hessian_and_gradient(const Vector &x, Matrix &H, Vector &g) const override {
            if (!unconstrained_->hessian_and_gradient(x, H, g)) {
                return false;
            }
            extend_hessian_and_gradient(x, H, g);
            return true;
        }

        void extend_hessian_and_gradient(const Vector &x, Matrix &H, Vector &g) const {
            Vector diff;

            if (box_->has_upper_bound()) {
                compute_diff_upper_bound(x, diff);
                g += current_barrier_parameter_ / diff;

                diff = pow2(diff);
                diff = current_barrier_parameter_ / diff;

                H.shift_diag(diff);
            }

            if (box_->has_lower_bound()) {
                compute_diff_lower_bound(x, diff);

                g += current_barrier_parameter_ / diff;

                diff = pow2(diff);
                diff = current_barrier_parameter_ / diff;

                H.shift_diag(diff);
            }
        }

        void compute_diff_upper_bound(const Vector &x, Vector &diff) const {
            if (diff.empty()) {
                diff.zeros(layout(x));
            }

            auto ub_view = local_view_device(*box_->upper_bound());
            auto x_view = local_view_device(x);
            auto diff_view = local_view_device(diff);

            auto zero = zero_;
            parallel_for(
                local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                    auto xi = x_view.get(i);
                    auto ubi = ub_view.get(i);
                    auto d = ubi - xi;
                    if (d == 0.) {
                        d = zero;
                    }

                    diff_view.set(i, d);
                });
        }

        void compute_diff_lower_bound(const Vector &x, Vector &diff) const {
            if (diff.empty()) {
                diff.zeros(layout(x));
            }

            auto lb_view = local_view_device(*box_->lower_bound());
            auto x_view = local_view_device(x);
            auto diff_view = local_view_device(diff);

            auto zero = zero_;
            parallel_for(
                local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                    auto xi = x_view.get(i);
                    auto lbi = lb_view.get(i);
                    auto d = xi - lbi;
                    if (d == 0.) {
                        d = zero;
                    }

                    diff_view.set(i, d);
                });
        }

        void extend_hessian(const Vector &x, Matrix &H) const {
            Vector diff;

            if (box_->has_upper_bound()) {
                compute_diff_upper_bound(x, diff);
                diff = pow2(diff);
                diff = current_barrier_parameter_ / diff;

                H.shift_diag(diff);
            }

            if (box_->has_lower_bound()) {
                compute_diff_lower_bound(x, diff);
                diff = pow2(diff);
                diff = current_barrier_parameter_ / diff;

                H.shift_diag(diff);
            }
        }

        void extend_gradient(const Vector &x, Vector &g) const {
            Vector diff;
            if (box_->has_upper_bound()) {
                compute_diff_upper_bound(x, diff);
                g += current_barrier_parameter_ / diff;
            }

            if (box_->has_lower_bound()) {
                compute_diff_lower_bound(x, diff);
                g -= current_barrier_parameter_ / diff;
            }
        }

        bool hessian_and_gradient(const Vector &x, Matrix &H, Matrix &preconditioner, Vector &g) const override {
            if (!unconstrained_->hessian_and_gradient(x, H, preconditioner, g)) {
                return false;
            }

            extend_hessian_and_gradient(x, H, g);
            return true;
        }

        bool has_preconditioner() const override { return unconstrained_->has_preconditioner(); }

        bool initialize_hessian(Matrix &H, Matrix &preconditioner) const override {
            if (!unconstrained_->initialize_hessian(H, preconditioner)) {
                return false;
            }

            return true;
        }

        bool value(const Vector &x, Scalar &value) const override {
            if (!unconstrained_->value(x, value)) {
                // return false;
                value = 0.0;
            }

            Scalar ub_value = 0.0;
            if (box_->has_upper_bound()) {
                ub_value = current_barrier_parameter_ * sum(logn(*box_->upper_bound() - x));
            }

            Scalar lb_value = 0.0;
            if (box_->has_lower_bound()) {
                ub_value = current_barrier_parameter_ * sum(logn(x - *box_->lower_bound()));
            }

            value -= (ub_value - lb_value);
            return true;
        }

        bool gradient(const Vector &x, Vector &g) const override {
            if (!unconstrained_->gradient(x, g)) {
                return false;
            }

            extend_gradient(x, g);
            return true;
        }

        bool update(const Vector &x) override {
            if (!unconstrained_->update(x)) {
                return false;
            }

            update_barrier();
            return true;
        }

        void set_barrier_parameter(const Scalar value) {
            barrier_parameter_ = value;
            current_barrier_parameter_ = value;
        }

        bool project_onto_feasibile_region(Vector &x) const override {
            // bool verbose = verbose_;
            if (box_->has_upper_bound()) {
                auto ub_view = local_view_device(*box_->upper_bound());
                auto x_view = local_view_device(x);

                Scalar soft_boundary = soft_boundary_;
                parallel_for(
                    local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                        auto xi = x_view.get(i);
                        auto ubi = ub_view.get(i);

                        if (xi > ubi) {
                            x_view.set(i, ubi - soft_boundary);
                        }
                    });
            }

            if (box_->has_lower_bound()) {
                auto lb_view = local_view_device(*box_->lower_bound());
                auto x_view = local_view_device(x);

                Scalar soft_boundary = soft_boundary_;
                parallel_for(
                    local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                        auto xi = x_view.get(i);
                        auto lbi = lb_view.get(i);

                        if (xi < lbi) {
                            x_view.set(i, lbi + soft_boundary);
                        }
                    });
            }

            return true;
        }

        void set_barrier_parameter_shrinking_factor(const Scalar value) { barrier_parameter_shrinking_factor_ = value; }

        void set_min_barrier_parameter(const Scalar value) { min_barrier_parameter_ = value; }

        void reset() { current_barrier_parameter_ = barrier_parameter_; }

    private:
        std::shared_ptr<Function> unconstrained_;
        std::shared_ptr<BoxConstraints> box_;

        Scalar barrier_parameter_{1e-10};
        Scalar barrier_parameter_shrinking_factor_{0.1};
        Scalar min_barrier_parameter_{1e-10};
        Scalar current_barrier_parameter_{1e-10};
        Scalar soft_boundary_{1e-7};
        Scalar zero_{1e-20};

        bool verbose_{false};

        void update_barrier() {
            current_barrier_parameter_ =
                std::max(current_barrier_parameter_ * barrier_parameter_shrinking_factor_, min_barrier_parameter_);

            if (verbose_) {
                utopia::out() << "current_barrier_parameter: " << current_barrier_parameter_ << '\n';
            }
        }
    };

}  // namespace utopia

#endif  // UTOPIA_LOG_BARRIER_FUNCTION_HPP
