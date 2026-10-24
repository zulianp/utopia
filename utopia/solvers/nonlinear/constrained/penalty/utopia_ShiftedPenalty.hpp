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
        using Layout = typename Traits<Vector>::Layout;
        using BoxConstraints = utopia::BoxConstraints<Vector>;
        using Transformation = utopia::Transformation<Matrix, Vector>;

        ShiftedPenalty() = default;
        explicit ShiftedPenalty(const std::shared_ptr<BoxConstraints> &box) : Super(box) {}
        virtual ~ShiftedPenalty() = default;

        std::string function_type() const override { return "ShiftedPenalty"; }

        // The penalty energy for the shift s and the active set a, both set by update(), is
        //   rho/2 sum_i dps_i(x)^2,  dps(x) = a (ub - x + s) (upper bound) or a (lb - x + s) (lower bound),
        // with gradient -rho dps(x) and hessian rho a. The active set is the one of the shift before the last
        // update. Determining it from the new shift instead makes Newton cycle between pushing the body out of the
        // constraint and releasing it
        void shifted_violation(const Vector &x, Vector &dps_out, Vector &active_out) const {
            Vector zero_shift;
            if (!shift) {
                zero_shift.zeros(layout(x));
            }

            const Vector &s = shift ? *shift : zero_shift;

            dps_out.zeros(layout(x));

            if (active) {
                active_out = *active;
            } else {
                active_out.zeros(layout(x));
            }

            auto s_view = const_local_view_device(s);
            auto x_view = const_local_view_device(x);
            auto dps_view = local_view_device(dps_out);
            auto a_view = local_view_device(active_out);

            // Without update() the active set is the one of the zero shift
            const bool has_active = static_cast<bool>(active);

            if (this->box()->has_upper_bound()) {
                auto ub_view = const_local_view_device(*this->box()->upper_bound());

                parallel_for(
                    local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                        const Scalar v = ub_view.get(i) - x_view.get(i) + s_view.get(i);
                        const Scalar a = has_active ? a_view.get(i) : Scalar(v <= 0);
                        dps_view.set(i, a * v);
                        a_view.set(i, a);
                    });
            } else if (this->box()->has_lower_bound()) {
                auto lb_view = const_local_view_device(*this->box()->lower_bound());

                parallel_for(
                    local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                        const Scalar v = lb_view.get(i) - x_view.get(i) + s_view.get(i);
                        const Scalar a = has_active ? a_view.get(i) : Scalar(v >= 0);
                        dps_view.set(i, a * v);
                        a_view.set(i, a);
                    });
            }
        }

        bool penalty_gradient(const Vector &x, Vector &g) const {
            const Scalar penalty_parameter = penalty_parameter_current_;

            Vector dps_x, active_x;
            shifted_violation(x, dps_x, active_x);

            auto dps_view = const_local_view_device(dps_x);
            auto g_view = local_view_device(g);

            parallel_for(
                local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                    g_view.set(i, g_view.get(i) - penalty_parameter * dps_view.get(i));
                });

            return true;
        }

        bool penalty_hessian(const Vector &x, Vector &H) const {
            const Scalar penalty_parameter = penalty_parameter_current_;

            if (H.empty()) {
                H.zeros(layout(x));
            }

            Vector dps_x, active_x;
            shifted_violation(x, dps_x, active_x);

            auto a_view = const_local_view_device(active_x);
            auto H_view = local_view_device(H);

            parallel_for(
                local_range_device(x), UTOPIA_LAMBDA(const SizeType i) { H_view.set(i, a_view.get(i) * penalty_parameter); });

            return true;
        }

        bool penalty_value(const Vector &x, Vector &v) const {
            const Scalar penalty_parameter = penalty_parameter_current_;

            Vector dps_x, active_x;
            shifted_violation(x, dps_x, active_x);

            auto dps_view = const_local_view_device(dps_x);
            auto v_view = local_view_device(v);

            parallel_for(
                local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                    const Scalar dps_i = dps_view.get(i);
                    v_view.set(i, 0.5 * penalty_parameter * dps_i * dps_i);
                });

            return true;
        }

        bool hessian_and_gradient(const Vector &x, Matrix &H, Vector &g) const override {
            return gradient(x, g) && hessian(x, H);
        }

        bool hessian(const Vector &x, Vector &H) const override {
            UTOPIA_TRACE_SCOPE("ShiftedPenalty::hessian");
            Vector work, h;

            if (this->has_transform()) {
                // Transform to constraint base
                this->transform()->transform(x, work);
                h.zeros(layout(work));
                if (!penalty_hessian(work, h)) {
                    return false;
                }

                this->apply_selection(h);

                // Scale contribuions
                if (this->has_scaling_matrix()) {
                    this->scaling_matrix()->apply(h, work);
                    h = work;
                }

                // The hessian in the constraint base is diag(h). In the problem base it is T^T diag(h) T,
                // whose diagonal is not T^T h
                Matrix diag_matrix = diag(h);
                Matrix temp;
                this->transform()->transform(diag_matrix, temp);
                Vector temp_diag = diag(temp);

                if (empty(H)) {
                    H = temp_diag;
                } else {
                    H += temp_diag;
                }
            } else {
                h.zeros(layout(x));
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
            Vector work, h;

            if (this->has_transform()) {
                // Transform to constraint base
                this->transform()->transform(x, work);
                h.zeros(layout(work));

                if (!penalty_hessian(work, h)) {
                    return false;
                }

                this->apply_selection(h);

                // Scale contribuions
                if (this->has_scaling_matrix()) {
                    this->scaling_matrix()->apply(h, work);
                    h = work;
                }

                work.zeros(layout(x));

                // Transform to problem base
                // this->transform()->inverse_transform_direction(h, work);
                // H.shift_diag(work);

                Matrix diag_matrix = diag(h);
                Matrix temp;
                this->transform()->transform(diag_matrix, temp);

                if (H.empty()) {
                    H = temp;
                } else {
                    H += temp;
                }

            } else {
                h.zeros(layout(x));
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

        // Adds the penalty energy to value
        bool value(const Vector &x, Scalar &value) const override {
            UTOPIA_TRACE_SCOPE("ShiftedPenalty::value");

            Vector work, v;

            if (this->has_transform()) {
                // Transform to constraint base
                this->transform()->transform(x, work);
                v.zeros(layout(work));

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
                value += sum(v);
            } else {
                v.zeros(layout(x));

                if (!penalty_value(x, v)) {
                    return false;
                }

                this->apply_selection(v);

                if (this->has_scaling_matrix()) {
                    this->scaling_matrix()->apply(v, work);
                    value += sum(work);
                } else {
                    value += sum(v);
                }
            }

            return true;
        }

        bool gradient(const Vector &x, Vector &grad) const override {
            UTOPIA_TRACE_SCOPE("ShiftedPenalty::gradient");
            Vector work, g;

            if (this->has_transform()) {
                // Transform to constraint base
                this->transform()->transform(x, work);
                g.zeros(layout(work));

                if (!penalty_gradient(work, g)) {
                    return false;
                }

                this->apply_selection(g);

                // Scale contribuions
                if (this->has_scaling_matrix()) {
                    this->scaling_matrix()->apply(g, work);
                    g = work;
                }

                // if (!layout(work).same(layout(grad))) {
                work.zeros(layout(grad));
                // }

                // Transform to problem base
                this->transform()->inverse_transform_direction(g, work);
                grad += work;
            } else {
                g.zeros(layout(x));
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

        SizeType count_active(const Vector &x) const {
            SizeType c = 0;
            if (this->box()->has_upper_bound()) {
                {
                    auto ub_view = const_local_view_device(*this->box()->upper_bound());
                    auto x_view = const_local_view_device(x);

                    parallel_reduce(
                        local_range_device(x),
                        UTOPIA_LAMBDA(const SizeType i)->Scalar { return x_view.get(i) > ub_view.get(i); },
                        c);
                }
            }

            c = x.comm().sum(c);
            return c;
        }

        // Augmented-Lagrangian update of the shift at x: s <- min(0, ub - x + s) (upper bound) or
        // s <- max(0, lb - x + s) (lower bound), and of the active set (where the shifted violation is not zero)
        void update_shift_aux(const Vector &x) {
            if (this->box()->has_upper_bound() && this->box()->has_lower_bound()) {
                if (!x.comm().rank()) {
                    m_utopia_warning_once(
                        "Upper and Lower bound simultaneously is not supported yet!\n"
                        "Enforcing upper bound only!\n");
                }
            }

            {
                auto x_view = const_local_view_device(x);
                auto s_view = local_view_device(*shift);
                auto a_view = local_view_device(*active);

                if (this->box()->has_upper_bound()) {
                    auto ub_view = const_local_view_device(*this->box()->upper_bound());

                    parallel_for(
                        local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                            const Scalar v = ub_view.get(i) - x_view.get(i) + s_view.get(i);
                            const bool a = v <= 0;
                            s_view.set(i, a * v);
                            a_view.set(i, a);
                        });
                } else if (this->box()->has_lower_bound()) {
                    auto lb_view = const_local_view_device(*this->box()->lower_bound());

                    parallel_for(
                        local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                            const Scalar v = lb_view.get(i) - x_view.get(i) + s_view.get(i);
                            const bool a = v >= 0;
                            s_view.set(i, a * v);
                            a_view.set(i, a);
                        });
                }
            }

            if (verbose_) {
                Vector active_x;
                shifted_violation(x, *dps, active_x);

                const SizeType c = sum(*active);
                Vector diff;
                if (this->box()->has_upper_bound()) {
                    diff = *this->box()->upper_bound() - x;
                    diff.e_min(0);
                } else {
                    diff = *this->box()->lower_bound() - x;
                    diff.e_max(0);
                }

                const Scalar penetration_norm = norm2(diff);
                const Scalar shift_norm = norm2(*shift);
                const Scalar dps_norm = norm2(*dps);

                if (!x.comm().rank() && c > 0) {
                    utopia::out() << (this->box()->has_upper_bound() ? "UB" : "LB") << ") Active: " << c
                                  << ", penetration_norm: " << penetration_norm << ", shift_norm: " << shift_norm
                                  << ", dps norm:" << dps_norm << "\n";
                }
            }
        }

        void reset_buffers(const Layout &l) {
            if (!shift) {
                active = std::make_shared<Vector>(l);
                shift = std::make_shared<Vector>(l);
                dps = std::make_shared<Vector>(l);
            } else if (!l.same(layout(*shift))) {
                active->zeros(l);
                shift->zeros(l);
                dps->zeros(l);
            }
        }

        void update_shift(const Vector &x) {
            UTOPIA_TRACE_SCOPE("ShiftedPenalty::shift");

            Vector work;
            if (this->has_transform()) {
                // Transform to constraint base
                this->transform()->transform(x, work);

                reset_buffers(layout(work));
                update_shift_aux(work);
            } else {
                reset_buffers(layout(x));
                update_shift_aux(x);
            }

            penalty_parameter_current_ =
                std::min(penalty_parameter_max_, penalty_parameter_factor_ * penalty_parameter_current_);
        }

        void update(const Vector &x) override { update_shift(x); }

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
            Super::read(in);

            penalty_parameter_max_ = 0;
            if (!Options()
                     .add_option(
                         "penalty_parameter", penalty_parameter_, "see Numerical Optimization - J. Nocedal, S. Wright.")
                     .add_option("penalty_parameter_max",
                                 penalty_parameter_max_,
                                 "Maximum value allowed for penalty_parameter.")
                     .add_option("penalty_parameter_factor",
                                 penalty_parameter_factor_,
                                 "Increase factor for the penalty parameter.")
                     .add_option("verbose",
                                 verbose_,
                                 "If verbose == true, prints additional information (with computational overhead).")
                     .add_option(
                         "zero_shift_on_reset", zero_shift_on_reset_, "Set the shift to zero when reset() is called.")
                     .parse(in)) {
                return;
            }

            penalty_parameter_current_ = penalty_parameter_;
            if (penalty_parameter_max_ == 0) {
                penalty_parameter_max_ = penalty_parameter_;
            }

            if (mpi_world_rank() == 0) {
                describe(utopia::out().stream());
            }
        }

        void reset() override {
            if (shift && zero_shift_on_reset_) {
                shift->set(0.);
            } else {
                if (shift && verbose_ && !zero_shift_on_reset_) {
                    Scalar norm_shift = norm2(*shift);
                    if (!shift->comm().rank()) {
                        utopia::out() << "ShiftedPenalty: keeping shift from previous work (mag: " << norm_shift
                                      << ")!\n";
                    }
                }
            }

            penalty_parameter_current_ = penalty_parameter_;
        }

        void describe(std::ostream &os) const {
            os << "-----------------------------------------\n";
            os << "utopia::ShiftedPenalty\n";
            os << "-----------------------------------------\n";
            os << "penalty_parameter:\t" << penalty_parameter_ << "\n";

            os << "-----------------------------------------\n";
        }

        void set_auxiliary_forcing(const std::shared_ptr<Vector> &vec) override { auxiliary_forcing_ = vec; }
        bool supports_auxiliary_forcing() const override { return true; }

        UTOPIA_NVCC_PRIVATE
        Scalar penalty_parameter_{1e4};
        Scalar penalty_parameter_current_{1e4};
        Scalar penalty_parameter_max_{1e4};
        Scalar penalty_parameter_factor_{1.5};

        std::shared_ptr<Vector> auxiliary_forcing_;
        std::shared_ptr<Vector> active;
        std::shared_ptr<Vector> shift;
        std::shared_ptr<Vector> dps;
        bool verbose_{false};
        bool zero_shift_on_reset_{true};
    };
}  // namespace utopia

#endif  // UTOPIA_SHIFTED_PENALTY_HPP
