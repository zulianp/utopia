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

        inline std::string function_type() const override { return "LogBarrierFunctionWithSelection"; }

        void read(Input &in) override {
            Super::read(in);
            in.get("infinity", infinity_);
            in.get("auto_selector", auto_selector_);
            in.get("skip_projection", skip_projection_);
            in.get("print_norms", print_norms_);
        }

        void ensure_selector() const {
            assert(boolean_selector_);
            if (boolean_selector_->empty()) {
                determine_boolean_selector();
            }
        }

        void allocate_selector() {
            if (!boolean_selector_) {
                boolean_selector_ = std::make_shared<Vector>();
            }
        }

        void determine_boolean_selector() const {
            if (!this->box_) return;

            if (this->verbose()) {
                utopia::out() << "LogBarrierFunctionWithSelection::determine_boolean_selector()\n";
            }

            if (empty(*boolean_selector_)) {
                if (this->box_->has_upper_bound()) {
                    boolean_selector_->zeros(layout(*this->box_->upper_bound()));
                } else if (this->box_->has_lower_bound()) {
                    boolean_selector_->zeros(layout(*this->box_->lower_bound()));
                }
            }

            if (this->box_->has_upper_bound() && this->box_->has_lower_bound()) {
                auto lb_view = local_view_device(*this->box_->lower_bound());
                auto ub_view = local_view_device(*this->box_->upper_bound());
                auto selector_view = local_view_device(*boolean_selector_);

                Scalar infinity = infinity_;
                parallel_for(
                    local_range_device(*boolean_selector_), UTOPIA_LAMBDA(const SizeType i) {
                        auto ubi = ub_view.get(i);
                        auto lbi = lb_view.get(i);

                        if (ubi < infinity || lbi > -infinity) {
                            selector_view.set(i, 1.);
                        }
                    });

            } else if (this->box_->has_upper_bound()) {
                auto ub_view = local_view_device(*this->box_->upper_bound());
                auto selector_view = local_view_device(*boolean_selector_);

                Scalar infinity = infinity_;
                parallel_for(
                    local_range_device(*boolean_selector_), UTOPIA_LAMBDA(const SizeType i) {
                        auto bound = ub_view.get(i);

                        if (bound < infinity) {
                            selector_view.set(i, 1.);
                        }
                    });
            } else if (this->box_->has_lower_bound()) {
                auto lb_view = local_view_device(*this->box_->lower_bound());
                auto selector_view = local_view_device(*boolean_selector_);

                Scalar infinity = infinity_;
                parallel_for(
                    local_range_device(*boolean_selector_), UTOPIA_LAMBDA(const SizeType i) {
                        auto bound = lb_view.get(i);

                        if (bound > -infinity) {
                            selector_view.set(i, 1.);
                        }
                    });
            }

            if (this->verbose()) {
                Scalar sum_selected = sum(*boolean_selector_);
                if (boolean_selector_->comm().rank() == 0) {
                    utopia::out() << "Selected: " << SizeType(sum_selected) << "\n";
                }
            }
        }

        void extend_hessian_and_gradient(const Vector &x, Matrix &H, Vector &g) const override {
            ensure_selector();

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
                g -= diff_selector;

                diff = pow2(diff);
                diff_selector = this->current_barrier_parameter_ / diff;
                diff_selector = e_mul(*boolean_selector_, diff);

                H.shift_diag(diff_selector);
            }
        }

        void extend_hessian(const Vector &x, Matrix &H) const override {
            ensure_selector();

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
            if (print_norms_) {
                Scalar norm_g = norm2(g);
                x.comm().root_print("Norm g");
                x.comm().root_print(norm_g);
            }

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

            if (print_norms_) {
                Scalar norm_g = norm2(g);
                x.comm().root_print("Norm g (with barrier)");
                x.comm().root_print(norm_g);
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

        bool extend_project_onto_feasibile_region(Vector &x) const override {
            if (skip_projection_) return true;

            assert(boolean_selector_);

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

            return true;
        }

        void set_selection(const std::shared_ptr<Vector> &boolean_selector) { boolean_selector_ = boolean_selector; }

        void reset() override {
            Super::reset();
            if (auto_selector_) {
                allocate_selector();
                determine_boolean_selector();
            }
        }

        void auto_selector(const bool value) { auto_selector_ = value; }

    private:
        bool auto_selector_{true};
        bool skip_projection_{false};
        bool print_norms_{false};
        std::shared_ptr<Vector> boolean_selector_;
        Scalar infinity_{std::numeric_limits<Scalar>::max()};
    };

}  // namespace utopia

#endif  // UTOPIA_LOG_BARRIER_FUNCTION_WITH_SELECTION_HPP
