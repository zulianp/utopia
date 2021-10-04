#ifndef UTOPIA_BOUNDED_LOG_BARRIER_FUNCTION_HPP
#define UTOPIA_BOUNDED_LOG_BARRIER_FUNCTION_HPP

#include "utopia_BoxConstraints.hpp"
#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_Layout.hpp"
#include "utopia_Options.hpp"

#include "utopia_LogBarrierFunctionBase.hpp"

#include <limits>

namespace utopia {
    template <class Matrix, class Vector>
    class BoundedLogBarrierFunction : public LogBarrierFunctionBase<Matrix, Vector> {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Function = utopia::Function<Matrix, Vector>;
        using BoxConstraints = utopia::BoxConstraints<Vector>;
        using Super = utopia::LogBarrierFunctionBase<Matrix, Vector>;

        BoundedLogBarrierFunction() {}

        BoundedLogBarrierFunction(const std::shared_ptr<Function> &unconstrained,
                                  const std::shared_ptr<BoxConstraints> &box)
            : Super(unconstrained, box) {}

        inline std::string function_type() const override { return "BoundedLogBarrierFunction"; }

        void read(Input &in) override {
            Super::read(in);
            if (!Options()
                     .add_option("barrier_thickness",
                                 barrier_thickness_,
                                 "see: Technical Supplement to Incremental Potential Contact: Intersection- and "
                                 "Inversion-free, Large-Deformation Dynamics.")
                     .parse(in)) {
                return;
            }
        }

        void add_barrier_gradient(const Vector &diff, Vector &g) const {
            auto diff_view = local_view_device(diff);
            auto g_view = local_view_device(g);

            auto d_hat = barrier_thickness_;

            // Currently it is not adaptive like in the paper
            auto stiffness = this->current_barrier_parameter_;

            parallel_for(
                local_range_device(diff), UTOPIA_LAMBDA(const SizeType i) {
                    auto d_i = diff_view.get(i);
                    auto d_m_d_hat = d_i - d_hat;

                    if (d_m_d_hat > 0) {
                        // Inside the thickness of the barrier
                        auto d_div_d_hat = d_i / d_hat;
                        auto b_g = 2 * device::log(d_div_d_hat) + (d_m_d_hat)*d_div_d_hat;
                        b_g *= -stiffness * d_m_d_hat;

                        auto g_i = g_view.get(i);
                        g_view.set(i, g_i + b_g);
                    }
                });
        }

        // !!! diff is modified inside !!!
        void add_barrier_hessian(Vector &diff, Matrix &hessian) const {
            {
                auto diff_view = local_view_device(diff);
                auto d_hat = barrier_thickness_;

                // Currently it is not adaptive like in the paper
                auto stiffness = this->current_barrier_parameter_;

                parallel_for(
                    local_range_device(diff), UTOPIA_LAMBDA(const SizeType i) {
                        auto d_i = diff_view.get(i);
                        auto d_m_d_hat = d_i - d_hat;

                        if (d_m_d_hat > 0) {
                            // Inside the thickness of the barrier
                            auto d_div_d_hat = d_i / d_hat;
                            auto a = 2 * device::log(d_div_d_hat) + (d_m_d_hat)*d_div_d_hat;
                            auto b = d_m_d_hat * (2 * d_div_d_hat + d_div_d_hat - d_m_d_hat * d_hat / (d_i * d_i));
                            auto b_H = -stiffness * (a + b);

                            diff_view.set(i, b_H);
                        }
                    });
            }

            hessian.shift_diag(diff);
        }

        void extend_hessian_and_gradient(const Vector &x, Matrix &H, Vector &g) const override {
            Vector diff;

            if (this->box_->has_upper_bound()) {
                this->compute_diff_upper_bound(x, diff);
                add_barrier_gradient(diff, g);
                add_barrier_hessian(diff, H);
            }

            if (this->box_->has_lower_bound()) {
                this->compute_diff_lower_bound(x, diff);
                add_barrier_gradient(diff, g);
                add_barrier_hessian(diff, H);
            }
        }

        void extend_hessian(const Vector &x, Matrix &H) const override {
            Vector diff;

            if (this->box_->has_upper_bound()) {
                this->compute_diff_upper_bound(x, diff);
                add_barrier_hessian(diff, H);
            }

            if (this->box_->has_lower_bound()) {
                this->compute_diff_lower_bound(x, diff);
                add_barrier_hessian(diff, H);
            }
        }

        void extend_gradient(const Vector &x, Vector &g) const override {
            Vector diff;

            if (this->box_->has_upper_bound()) {
                this->compute_diff_upper_bound(x, diff);
                add_barrier_gradient(diff, g);
            }

            if (this->box_->has_lower_bound()) {
                this->compute_diff_lower_bound(x, diff);
                add_barrier_gradient(diff, g);
            }
        }

        void add_barrier_value(const Vector &diff, Scalar &val) const {
            auto diff_view = local_view_device(diff);

            auto d_hat = barrier_thickness_;

            // Currently it is not adaptive like in the paper
            auto stiffness = this->current_barrier_parameter_;

            Scalar b_val = 0.;
            parallel_reduce(
                local_range_device(diff),
                UTOPIA_LAMBDA(const SizeType i) {
                    auto d_i = diff_view.get(i);
                    auto d_m_d_hat = d_i - d_hat;

                    if (d_m_d_hat > 0) {
                        auto d_div_d_hat = d_i / d_hat;
                        // Inside the thickness of the barrier
                        return -stiffness * d_m_d_hat * d_m_d_hat * device::log(d_div_d_hat);
                    } else {
                        return 0.0;
                    }
                },
                b_val);

            val += b_val;
        }

        void extend_value(const Vector &x, Scalar &value) const override {
            bool must_all_reduce = false;

            Vector diff;

            Scalar ub_value = 0.0;
            if (this->box_->has_upper_bound()) {
                must_all_reduce = true;
                diff = *this->box_->upper_bound() - x;
                add_barrier_value(diff, ub_value);
            }

            Scalar lb_value = 0.0;
            if (this->box_->has_lower_bound()) {
                must_all_reduce = true;
                diff = x - *this->box_->lower_bound();
                add_barrier_value(diff, lb_value);
            }

            if (must_all_reduce) {
                auto global_value = x.comm().sum(lb_value + ub_value);
                value += global_value;
            }
        }

        bool project_onto_feasibile_region(Vector &x) const override {
            // bool verbose = verbose_;
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

    private:
        Scalar barrier_thickness_{1e-4};
    };

}  // namespace utopia

#endif  // UTOPIA_BOUNDED_LOG_BARRIER_FUNCTION_HPP
