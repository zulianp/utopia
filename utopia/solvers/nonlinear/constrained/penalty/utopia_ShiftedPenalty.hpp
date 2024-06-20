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

        bool penalty_gradient(const Vector &x, Vector &g) const {
            const Scalar penalty_parameter = penalty_parameter_current_;

            auto dps_view = const_local_view_device(*dps);
            auto g_view = local_view_device(g);

            // SizeType cc = count_active(x);
            // if (!x.comm().rank()) {
            //     utopia::out() << "Actual active " << cc << "\n";
            // }

            parallel_for(
                local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                    const Scalar dps = dps_view.get(i);
                    g_view.set(i, g_view.get(i) - penalty_parameter * dps);
                });

            return true;
        }

        bool penalty_hessian(const Vector &x, Vector &H) const {
            Vector d, diag_B;
            const Scalar penalty_parameter = penalty_parameter_current_;

            if (H.empty()) {
                H.zeros(layout(x));
            }

            if (!active) {
                Utopia::Abort();
            }

            auto a_view = const_local_view_device(*active);
            auto diag_B_view = local_view_device(H);

            parallel_for(
                local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                    const Scalar ai = a_view.get(i);
                    diag_B_view.set(i, ai * penalty_parameter);
                });

            return true;
        }

        bool penalty_value(const Vector &x, Vector &v) const {
            Vector d;
            const Scalar penalty_parameter = penalty_parameter_current_;

            auto a_view = local_view_device(*active);

            if (this->box()->has_upper_bound()) {
                d = *this->box()->upper_bound() - x;

                auto d_view = const_local_view_device(d);
                auto v_view = local_view_device(v);

                parallel_for(
                    local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                        const Scalar di = d_view.get(i);
                        const Scalar active = a_view.get(i);
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
                        const Scalar active = a_view.get(i);
                        v_view.set(i, active * di * di * penalty_parameter);
                    });
            }

            // if (this->box()->has_upper_bound() || this->box()->has_lower_bound()) {
            //     auto d_view = const_local_view_device(*dps);
            //     auto v_view = local_view_device(v);

            //     parallel_for(
            //         local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
            //             const Scalar di = d_view.get(i);
            //             v_view.set(i, di * di * penalty_parameter);
            //         });
            // }

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

                work.zeros(layout(H));

                // Transform to problem base
                this->transform()->inverse_transform_direction(h, work);
                H += (work);
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

        /// FIXME shift is not considred in the energy value here
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

        void update_shift_aux(const Vector &x) {
            if (this->box()->has_upper_bound() && this->box()->has_lower_bound()) {
                if (!x.comm().rank()) {
                    m_utopia_warning_once(
                        "Upper and Lower bound simultaneously is not supported yet!\n"
                        "Enforcing upper bound only!\n");
                }
            }

            if (this->box()->has_upper_bound()) {
                {
                    // {
                    //     std::stringstream ss;
                    //     ss << "x: " << x.local_size() << " " << x.size() << "\n";
                    //     ss << "shift: " << shift->local_size() << " " << shift->size() << "\n";
                    //     ss << "active: " << active->local_size() << " " << active->size() << "\n";
                    //     ss << "dps: " << dps->local_size() << " " << dps->size() << "\n";
                    //     ss << "upper_bound: " << this->box()->upper_bound()->local_size() << " "
                    //        << this->box()->upper_bound()->size() << "\n";

                    //     x.comm().synched_print(ss.str());
                    // }

                    auto ub_view = const_local_view_device(*this->box()->upper_bound());
                    auto x_view = const_local_view_device(x);
                    auto s_view = local_view_device(*shift);
                    auto a_view = local_view_device(*active);
                    auto dps_view = local_view_device(*dps);

                    parallel_for(
                        local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                            const auto ubi = ub_view.get(i);
                            const auto xi = x_view.get(i);
                            const auto di = ubi - xi;
                            auto si = s_view.get(i);

                            auto dps = di + si;
                            const bool active = dps <= 0;

                            si = active * dps;
                            dps = active * (di + si);

                            a_view.set(i, active);
                            s_view.set(i, si);
                            dps_view.set(i, dps);

                            // printf("%d) %g %g\n", i, si, dps);
                        });
                }

                if (verbose_) {
                    const SizeType c = sum(*active);
                    Vector diff = *this->box()->upper_bound() - x;
                    diff.e_min(0);
                    const Scalar penetration_norm = norm2(diff);
                    const Scalar shift_norm = norm2(*shift);

                    if (!x.comm().rank() && c > 0) {
                        utopia::out() << "UB) Active: " << c << ", penetration_norm: " << penetration_norm
                                      << ", shift_norm: " << shift_norm << "\n";
                    }
                }
            } else if (this->box()->has_lower_bound()) {
                {
                    auto lb_view = const_local_view_device(*this->box()->lower_bound());
                    auto x_view = const_local_view_device(x);
                    auto s_view = local_view_device(*shift);
                    auto a_view = local_view_device(*active);
                    auto dps_view = local_view_device(*dps);

                    parallel_for(
                        local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                            const auto lbi = lb_view.get(i);
                            const auto xi = x_view.get(i);
                            const auto di = lbi - xi;
                            auto si = s_view.get(i);

                            auto dps = di + si;
                            const bool active = dps >= 0;

                            si = active * dps;
                            dps = active * (di + si);

                            a_view.set(i, active);
                            s_view.set(i, si);
                            dps_view.set(i, dps);
                        });
                }

                if (verbose_) {
                    const SizeType c = sum(*active);
                    Vector diff = *this->box()->lower_bound() - x;
                    diff.e_max(0);
                    const Scalar penetration_norm = norm2(diff);
                    const Scalar shift_norm = norm2(*shift);
                    const Scalar dps_norm = norm2(*dps);

                    if (!x.comm().rank() && c > 0) {
                        utopia::out() << "LB) Active: " << c << ", penetration_norm: " << penetration_norm
                                      << ", shift_norm: " << shift_norm << ", dps norm:" << dps_norm << "\n";
                    }
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
            if (shift) {
                shift->set(0.);
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
    };
}  // namespace utopia

#endif  // UTOPIA_SHIFTED_PENALTY_HPP
