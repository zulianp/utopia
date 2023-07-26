#ifndef UTOPIA_SHIFTED_PENALTY_HPP
#define UTOPIA_SHIFTED_PENALTY_HPP

#include "utopia_BoxConstraints.hpp"
#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_Layout.hpp"
#include "utopia_Options.hpp"

#include "utopia_Penalty.hpp"

#include <limits>

namespace utopia {

    template <class Matrix, class Vector = typename Traits<Matrix>::Vector>
    class ShiftedPenalty : public Penalty<Matrix, Vector> {
    public:
        using Super = utopia::Penalty<Matrix, Vector>;
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using BoxConstraints = utopia::BoxConstraints<Vector>;
        using Transformation = utopia::Transformation<Matrix, Vector>;

        ShiftedPenalty() = default;
        explicit ShiftedPenalty(const std::shared_ptr<BoxConstraints> &box) : Super(box) {}
        virtual ~ShiftedPenalty() = default;

        std::string function_type() const override { return "ShiftedPenalty"; }

        bool penalty_gradient(const Vector &x, Vector &g) const {
            Vector d;
            const Scalar penalty_parameter = penalty_parameter_;

            if (this->box()->has_upper_bound()) {
                d = *this->box()->upper_bound() - x;

                auto d_view = const_local_view_device(d);
                auto g_view = local_view_device(g);

                parallel_for(
                    local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                        const Scalar di = d_view.get(i);
                        const Scalar gi = g_view.get(i);
                        const Scalar active = di <= 0;
                        const Scalar g_active = active * (gi - penalty_parameter * di);
                        const Scalar gg = -(gi + g_active);
                        g_view.set(i, gg);
                        // diag_B_view.set(i, active * penalty_parameter);
                    });
            }

            if (this->box()->has_lower_bound()) {
                // TODO Check
                d = *this->box()->lower_bound() - x;

                auto d_view = const_local_view_device(d);
                auto g_view = local_view_device(g);

                parallel_for(
                    local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                        const Scalar di = d_view.get(i);
                        const Scalar gi = g_view.get(i);
                        const Scalar active = di >= 0;
                        const Scalar g_active = active * (gi - penalty_parameter * di);
                        const Scalar gg = -(gi + g_active);
                        g_view.set(i, gg);
                        // diag_B_view.set(i, active * penalty_parameter);
                    });
            }

            return true;
        }

        bool penalty_hessian(const Vector &x, Vector &H) const {
            Vector d, diag_B;
            const Scalar penalty_parameter = penalty_parameter_;

            if (H.empty()) {
                H.zeros(layout(x));
            }

            if (this->box()->has_upper_bound()) {
                d = *this->box()->upper_bound() - x;

                auto d_view = const_local_view_device(d);
                auto diag_B_view = local_view_device(H);

                parallel_for(
                    local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                        const Scalar di = d_view.get(i);
                        const Scalar active = di <= 0;

                        diag_B_view.set(i, active * penalty_parameter);
                    });
            }

            if (this->box()->has_lower_bound()) {
                // TODO Check
                d = *this->box()->lower_bound() - x;

                auto d_view = const_local_view_device(d);
                auto diag_B_view = local_view_device(H);

                parallel_for(
                    local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                        const Scalar di = d_view.get(i);
                        const Scalar active = di >= 0;

                        diag_B_view.set(i, active * penalty_parameter);
                    });
            }

            return true;
        }

        bool penalty_value(const Vector &x, Vector &v) const {
            Vector d;
            const Scalar penalty_parameter = penalty_parameter_;

            if (this->box()->has_upper_bound()) {
                d = *this->box()->upper_bound() - x;

                auto d_view = const_local_view_device(d);
                auto v_view = local_view_device(v);

                parallel_for(
                    local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                        const Scalar di = d_view.get(i);
                        const Scalar active = di <= 0;
                        v_view.set(i, active * di * di * penalty_parameter);
                    });
            }

            if (this->box()->has_lower_bound()) {
                // TODO Check
                d = *this->box()->lower_bound() - x;

                auto d_view = const_local_view_device(d);
                auto v_view = local_view_device(v);

                parallel_for(
                    local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                        const Scalar di = d_view.get(i);
                        const Scalar active = di >= 0;
                        v_view.set(i, active * di * di * penalty_parameter);
                    });
            }

            return true;
        }

        bool hessian_and_gradient(const Vector &x, Matrix &H, Vector &g) const override {
            return gradient(x, g) && hessian(x, H);
        }

        bool hessian(const Vector &x, Vector &H) const override {
            UTOPIA_TRACE_SCOPE("ShiftedPenalty::hessian");
            Vector work, h(layout(x), 0.);

            if (this->has_transform()) {
                // Transform to constraint base
                this->transform()->transform(x, work);
                if (!penalty_hessian(work, h)) {
                    return false;
                }

                this->apply_selection(h);

                // Scale contribuions
                if (this->has_scaling_matrix()) {
                    this->scaling_matrix()->apply(h, work);
                    h = work;
                }

                // Transform to problem base
                this->transform()->inverse_transform_direction(h, work);
                H += (work);
            } else {
                if (!penalty_hessian(x, h)) {
                    return false;
                }

                this->apply_selection(h);

                if (this->has_scaling_matrix()) {
                    this->scaling_matrix()->apply(h, work);
                    if (empty(H)) {
                        H = work;
                    } else {
                        H += (work);
                    }
                } else {
                    if (empty(H)) {
                        H = h;
                    } else {
                        H += (h);
                    }
                }
            }

            return true;
        }

        bool hessian(const Vector &x, Matrix &H) const override {
            UTOPIA_TRACE_SCOPE("ShiftedPenalty::hessian");
            Vector work, h(layout(x), 0.);

            if (this->has_transform()) {
                // Transform to constraint base
                this->transform()->transform(x, work);
                if (!penalty_hessian(work, h)) {
                    return false;
                }

                this->apply_selection(h);

                // Scale contribuions
                if (this->has_scaling_matrix()) {
                    this->scaling_matrix()->apply(h, work);
                    h = work;
                }

                // Transform to problem base
                this->transform()->inverse_transform_direction(h, work);
                H.shift_diag(work);
            } else {
                if (!penalty_hessian(x, h)) {
                    return false;
                }

                this->apply_selection(h);

                if (this->has_scaling_matrix()) {
                    this->scaling_matrix()->apply(h, work);
                    H.shift_diag(work);
                } else {
                    H.shift_diag(h);
                }
            }

            return true;
        }

        /// FIXME shift is not considred in the energy value here
        bool value(const Vector &x, Scalar &value) const override {
            UTOPIA_TRACE_SCOPE("ShiftedPenalty::value");

            Vector work, v(layout(x), 0.);

            if (this->has_transform()) {
                // Transform to constraint base
                this->transform()->transform(x, work);
                if (!penalty_value(work, v)) {
                    return false;
                }

                this->apply_selection(v);

                // Scale contribuions
                if (this->has_scaling_matrix()) {
                    this->scaling_matrix()->apply(v, work);
                    v = work;
                }

                // Transform to problem base
                this->transform()->inverse_transform_direction(v, work);
                value = sum(v);
            } else {
                if (!penalty_value(x, v)) {
                    return false;
                }

                this->apply_selection(v);

                if (this->has_scaling_matrix()) {
                    this->scaling_matrix()->apply(v, work);
                    value = sum(work);
                } else {
                    value = sum(v);
                }
            }

            return true;
        }

        bool gradient(const Vector &x, Vector &grad) const override {
            UTOPIA_TRACE_SCOPE("ShiftedPenalty::gradient");
            Vector work, g(layout(x), 0.);

            if (this->has_transform()) {
                // Transform to constraint base
                this->transform()->transform(x, work);

                // Transform gradient
                this->transform()->transform(grad, g);

                if (!penalty_gradient(work, g)) {
                    return false;
                }

                this->apply_selection(g);

                // Scale contribuions
                if (this->has_scaling_matrix()) {
                    this->scaling_matrix()->apply(g, work);
                    g = work;
                }

                // Transform to problem base
                this->transform()->inverse_transform_direction(g, work);
                grad += work;
            } else {
                // Copy gradient
                g = grad;
                if (!penalty_gradient(x, g)) {
                    return false;
                }

                this->apply_selection(g);

                if (this->has_scaling_matrix()) {
                    this->scaling_matrix()->apply(g, work);
                    grad += work;
                } else {
                    grad += g;
                }
            }

            return true;
        }

        void compute_diff_upper_bound(const Vector &x, Vector &diff) const {
            if (diff.empty()) {
                diff.zeros(layout(x));
            }

            auto ub_view = local_view_device(*this->box()->upper_bound());
            auto x_view = local_view_device(x);
            auto diff_view = local_view_device(diff);

            parallel_for(
                local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                    auto xi = x_view.get(i);
                    auto ubi = ub_view.get(i);
                    auto d = ubi - xi;
                    diff_view.set(i, d);
                });
        }

        void compute_diff_lower_bound(const Vector &x, Vector &diff) const {
            if (diff.empty()) {
                diff.zeros(layout(x));
            }

            auto lb_view = local_view_device(*this->box()->lower_bound());
            auto x_view = local_view_device(x);
            auto diff_view = local_view_device(diff);

            parallel_for(
                local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                    auto xi = x_view.get(i);
                    auto lbi = lb_view.get(i);
                    auto d = xi - lbi;
                    diff_view.set(i, d);
                });
        }

        void read(Input &in) override {
            if (!Options()
                     .add_option(
                         "penalty_parameter", penalty_parameter_, "see Numerical Optimization - J. Nocedal, S. Wright.")
                     .parse(in)) {
                return;
            }

            if (mpi_world_rank() == 0) {
                describe(utopia::out().stream());
            }
        }

        void describe(std::ostream &os) const {
            os << "-----------------------------------------\n";
            os << "utopia::ShiftedPenalty\n";
            os << "-----------------------------------------\n";
            os << "penalty_parameter:\t" << penalty_parameter_ << "\n";

            os << "-----------------------------------------\n";
        }

        UTOPIA_NVCC_PRIVATE
        Scalar penalty_parameter_{1e8};
    };
}  // namespace utopia

#endif  // UTOPIA_SHIFTED_PENALTY_HPP
