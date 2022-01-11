#ifndef UTOPIA_LOG_BARRIER_FUNCTION_HPP
#define UTOPIA_LOG_BARRIER_FUNCTION_HPP

#include "utopia_make_unique.hpp"

#include "utopia_BoxConstraints.hpp"
#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_Layout.hpp"
#include "utopia_Options.hpp"

#include "utopia_LogBarrierFunctionBase.hpp"

#include <limits>

namespace utopia {

    template <class Matrix, class Vector>
    class LogBarrier : public LogBarrierBase<Matrix, Vector> {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Function = utopia::Function<Matrix, Vector>;
        using BoxConstraints = utopia::BoxConstraints<Vector>;
        using Super = utopia::LogBarrierBase<Matrix, Vector>;

        LogBarrier() = default;
        explicit LogBarrier(const std::shared_ptr<BoxConstraints> &box) : Super(box) {}

        void hessian_and_gradient(const Vector &x, Matrix &H, Vector &g) const override {
            Vector diff;

            if (this->box_->has_upper_bound()) {
                this->compute_diff_upper_bound(x, diff);

                if (this->scaling_matrix()) {
                    g += (*this->scaling_matrix()) * (this->current_barrier_parameter_ / diff);
                } else {
                    g += this->current_barrier_parameter_ / diff;
                }

                diff = pow2(diff);
                diff = this->current_barrier_parameter_ / diff;

                if (this->scaling_matrix()) {
                    H.shift_diag((*this->scaling_matrix()) * diff);
                } else {
                    H.shift_diag(diff);
                }
            }

            if (this->box_->has_lower_bound()) {
                this->compute_diff_lower_bound(x, diff);

                if (this->scaling_matrix()) {
                    g += (*this->scaling_matrix()) * (this->current_barrier_parameter_ / diff);
                } else {
                    g += this->current_barrier_parameter_ / diff;
                }

                diff = pow2(diff);
                diff = this->current_barrier_parameter_ / diff;

                if (this->scaling_matrix()) {
                    H.shift_diag((*this->scaling_matrix()) * diff);
                } else {
                    H.shift_diag(diff);
                }
            }
        }

        void hessian(const Vector &x, Matrix &H) const override {
            Vector diff;

            if (this->box_->has_upper_bound()) {
                this->compute_diff_upper_bound(x, diff);
                diff = pow2(diff);
                diff = this->current_barrier_parameter_ / diff;

                if (this->scaling_matrix()) {
                    H.shift_diag((*this->scaling_matrix()) * diff);
                } else {
                    H.shift_diag(diff);
                }
            }

            if (this->box_->has_lower_bound()) {
                this->compute_diff_lower_bound(x, diff);
                diff = pow2(diff);
                diff = this->current_barrier_parameter_ / diff;

                if (this->scaling_matrix()) {
                    H.shift_diag((*this->scaling_matrix()) * diff);
                } else {
                    H.shift_diag(diff);
                }
            }
        }

        void hessian_diag(const Vector &x, Vector &h) const override {
            // TODO
            Utopia::Abort("IMPLEMENT ME");
        }

        void gradient(const Vector &x, Vector &g) const override {
            Vector diff;
            if (this->box_->has_upper_bound()) {
                this->compute_diff_upper_bound(x, diff);

                if (this->scaling_matrix()) {
                    g += (*this->scaling_matrix()) * (this->current_barrier_parameter_ / diff);
                } else {
                    g += this->current_barrier_parameter_ / diff;
                }
            }

            if (this->box_->has_lower_bound()) {
                this->compute_diff_lower_bound(x, diff);

                if (this->scaling_matrix()) {
                    g -= (*this->scaling_matrix()) * (this->current_barrier_parameter_ / diff);
                } else {
                    g -= this->current_barrier_parameter_ / diff;
                }
            }
        }

        void value(const Vector &x, Scalar &value) const override {
            Scalar ub_value = 0.0;
            if (this->box_->has_upper_bound()) {
                if (this->scaling_matrix()) {
                    ub_value = this->current_barrier_parameter_ *
                               sum((*this->scaling_matrix()) * logn(*this->box_->upper_bound() - x));
                } else {
                    ub_value = this->current_barrier_parameter_ * sum(logn(*this->box_->upper_bound() - x));
                }
            }

            Scalar lb_value = 0.0;
            if (this->box_->has_lower_bound()) {
                if (this->scaling_matrix()) {
                    ub_value = this->current_barrier_parameter_ *
                               sum((*this->scaling_matrix()) * logn(x - *this->box_->lower_bound()));
                } else {
                    ub_value = this->current_barrier_parameter_ * sum(logn(x - *this->box_->lower_bound()));
                }
            }

            value -= (ub_value - lb_value);
        }

        bool project_onto_feasibile_region(Vector &x) const override {
            if (this->box_->has_upper_bound()) {
                auto ub_view = local_view_device(*this->box_->upper_bound());
                auto x_view = local_view_device(x);

                Scalar soft_boundary = this->soft_boundary_;
                parallel_for(
                    local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                        auto xi = x_view.get(i);
                        auto ubi = ub_view.get(i);

                        if (xi > ubi) {
                            x_view.set(i, ubi - soft_boundary);
                        }
                    });
            }

            if (this->box_->has_lower_bound()) {
                auto lb_view = local_view_device(*this->box_->lower_bound());
                auto x_view = local_view_device(x);

                Scalar soft_boundary = this->soft_boundary_;
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

        void read(Input &in) override { Super::read(in); }
    };

    template <class Matrix, class Vector>
    class LogBarrierFunction : public LogBarrierFunctionBase<Matrix, Vector> {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Function = utopia::Function<Matrix, Vector>;
        using BoxConstraints = utopia::BoxConstraints<Vector>;
        using Super = utopia::LogBarrierFunctionBase<Matrix, Vector>;
        using LogBarrier = utopia::LogBarrier<Matrix, Vector>;

        LogBarrierFunction() { this->set_barrier(std::make_shared<LogBarrier>()); }

        LogBarrierFunction(const std::shared_ptr<Function> &unconstrained, const std::shared_ptr<BoxConstraints> &box)
            : Super(unconstrained, std::make_shared<LogBarrier>(box)) {}

        inline std::string function_type() const override { return "LogBarrierFunction"; }
    };

}  // namespace utopia

#endif  // UTOPIA_LOG_BARRIER_FUNCTION_HPP
