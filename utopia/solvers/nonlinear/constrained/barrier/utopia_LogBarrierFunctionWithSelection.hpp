#ifndef UTOPIA_LOG_BARRIER_FUNCTION_WITH_SELECTION_HPP
#define UTOPIA_LOG_BARRIER_FUNCTION_WITH_SELECTION_HPP

#include "utopia_BoxConstraints.hpp"
#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_Layout.hpp"
#include "utopia_Options.hpp"

#include "utopia_LogBarrierFunctionBase.hpp"

#include <limits>

namespace utopia {
    template <class Matrix, class Vector>
    class LogBarrierFunctionWithSelection : public LogBarrierFunctionBase<Matrix, Vector> {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Function = utopia::Function<Matrix, Vector>;
        using BoxConstraints = utopia::BoxConstraints<Vector>;
        using Super = utopia::LogBarrierFunctionBase<Matrix, Vector>;

        LogBarrierFunctionWithSelection() {}

        LogBarrierFunctionWithSelection(const std::shared_ptr<Function> &unconstrained,
                                        const std::shared_ptr<BoxConstraints> &box)
            : Super(unconstrained, box) {}

        void extend_hessian_and_gradient(const Vector &x, Matrix &H, Vector &g) const override {
            Vector diff, diff_selector;

            if (this->box_->has_upper_bound()) {
                this->compute_diff_upper_bound(x, diff);
                diff_selector = this->current_barrier_parameter_ / diff;
                diff_selector = e_mul(*boolean_selector_, diff);
                g += diff_selector;

                diff = pow2(diff);
                diff_selector = this->current_barrier_parameter_ / diff;
                diff_selector = e_mul(*boolean_selector_, diff);

                H.shift_diag(diff_selector);
            }

            if (this->box_->has_lower_bound()) {
                this->compute_diff_lower_bound(x, diff);

                diff_selector = this->current_barrier_parameter_ / diff;
                diff_selector = e_mul(*boolean_selector_, diff);
                g += diff_selector;

                diff = pow2(diff);
                diff_selector = this->current_barrier_parameter_ / diff;
                diff_selector = e_mul(*boolean_selector_, diff);

                H.shift_diag(diff_selector);
            }
        }

        void extend_hessian(const Vector &x, Matrix &H) const override {
            Vector diff;

            if (this->box_->has_upper_bound()) {
                this->compute_diff_upper_bound(x, diff);
                diff = pow2(diff);
                diff = this->current_barrier_parameter_ / diff;
                diff = e_mul(*boolean_selector_, diff);

                H.shift_diag(diff);
            }

            if (this->box_->has_lower_bound()) {
                this->compute_diff_lower_bound(x, diff);
                diff = pow2(diff);
                diff = this->current_barrier_parameter_ / diff;
                diff = e_mul(*boolean_selector_, diff);

                H.shift_diag(diff);
            }
        }

        void extend_gradient(const Vector &x, Vector &g) const override {
            Vector diff;
            if (this->box_->has_upper_bound()) {
                this->compute_diff_upper_bound(x, diff);
                diff = this->current_barrier_parameter_ / diff;
                diff = e_mul(*boolean_selector_, diff);
                g += diff;
            }

            if (this->box_->has_lower_bound()) {
                this->compute_diff_lower_bound(x, diff);
                diff = this->current_barrier_parameter_ / diff;
                diff = e_mul(*boolean_selector_, diff);
                g -= diff;
            }
        }

        void extend_value(const Vector &x, Scalar &value) const override {
            Scalar ub_value = 0.0;
            if (this->box_->has_upper_bound()) {
                ub_value = this->current_barrier_parameter_ *
                           sum(e_mul(*boolean_selector_, logn(*this->box_->upper_bound() - x)));
            }

            Scalar lb_value = 0.0;
            if (this->box_->has_lower_bound()) {
                ub_value = this->current_barrier_parameter_ *
                           sum(e_mul(*boolean_selector_, logn(x - *this->box_->lower_bound())));
            }

            value -= (ub_value - lb_value);
        }

        bool project_onto_feasibile_region(Vector &x) const override {
            // bool verbose = verbose_;
            if (this->box_->has_upper_bound()) {
                auto ub_view = local_view_device(*this->box_->upper_bound());
                auto x_view = local_view_device(x);
                auto selector_view = local_view_device(*boolean_selector_);

                Scalar soft_boundary = this->soft_boundary_;
                parallel_for(
                    local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                        auto s = selector_view.get(i);

                        if (s > 0.99) {
                            auto xi = x_view.get(i);
                            auto ubi = ub_view.get(i);

                            if (xi > ubi) {
                                x_view.set(i, ubi - soft_boundary);
                            }
                        }
                    });
            }

            if (this->box_->has_lower_bound()) {
                auto lb_view = local_view_device(*this->box_->lower_bound());
                auto x_view = local_view_device(x);
                auto selector_view = local_view_device(*boolean_selector_);

                Scalar soft_boundary = this->soft_boundary_;
                parallel_for(
                    local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                        auto s = selector_view.get(i);

                        if (s > 0.99) {
                            auto xi = x_view.get(i);
                            auto lbi = lb_view.get(i);

                            if (xi < lbi) {
                                x_view.set(i, lbi + soft_boundary);
                            }
                        }
                    });
            }
        }

        void set_selector(const std::shared_ptr<Vector> &boolean_selector) { boolean_selector_ = boolean_selector; }

    private:
        std::shared_ptr<Vector> boolean_selector_;
    };

}  // namespace utopia

#endif  // UTOPIA_LOG_BARRIER_FUNCTION_WITH_SELECTION_HPP
