#ifndef UTOPIA_SHIFTED_PENALTY_HPP
#define UTOPIA_SHIFTED_PENALTY_HPP

#include "utopia_BoxConstraints.hpp"
#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_Layout.hpp"
#include "utopia_Options.hpp"

#include <limits>

namespace utopia {

    template <class Matrix, class Vector = typename Traits<Matrix>::Vector>
    class ShiftedPenalty : public Configurable {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using BoxConstraints = utopia::BoxConstraints<Vector>;
        using Transformation = utopia::Transformation<Matrix, Vector>;

        ShiftedPenalty() = default;
        explicit ShiftedPenalty(const std::shared_ptr<BoxConstraints> &box) : box_(box) {}
        virtual ~ShiftedPenalty() = default;

        bool penalty_gradient(const Vector &x, Vector &g) const {
            Vector d;
            const Scalar penalty_parameter = penalty_parameter_;

            if (box_->has_upper_bound()) {
                d = *box_->upper_bound() - x;

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

            if (box_->has_lower_bound()) {
                // TODO Check
                d = *box_->lower_bound() - x;

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
        }

        bool penalty_hessian(const Vector &x, Vector &H) const {
            Vector d, diag_B;
            const Scalar penalty_parameter = penalty_parameter_;

            diag_B.zeros(layout(x));

            if (box_->has_upper_bound()) {
                d = *box_->upper_bound() - x;

                auto d_view = const_local_view_device(d);
                auto diag_B_view = local_view_device(diag_B);

                parallel_for(
                    local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                        const Scalar di = d_view.get(i);
                        const Scalar active = di <= 0;

                        diag_B_view.set(i, active * penalty_parameter);
                    });
            }

            if (box_->has_lower_bound()) {
                // TODO Check
                d = *box_->lower_bound() - x;

                auto d_view = const_local_view_device(d);
                auto diag_B_view = local_view_device(diag_B);

                parallel_for(
                    local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                        const Scalar di = d_view.get(i);
                        const Scalar active = di >= 0;

                        diag_B_view.set(i, active * penalty_parameter);
                    });
            }

            H.shift_diag(diag_B);
        }

        bool penalty_value(const Vector &x, Vector &v) const {
            const Scalar penalty_parameter = penalty_parameter_;
            //
        }

        bool hessian_and_gradient(const Vector &x, Matrix &H, Vector &g) const {
            return gradient(x, g) && hessian(x, H);
        }

        bool hessian(const Vector &x, Vector &H) const {
            Vector work, h(layout(x), 0.);

            if (has_transform()) {
                // Transform to constraint base
                transform_->transform(x, work);
                if (!penalty_hessian(work, h)) {
                    return false;
                }

                apply_selection(h);

                // Scale contribuions
                if (has_scaling_matrix()) {
                    scaling_matrix()->apply(h, work);
                    h = work;
                }

                // Transform to problem base
                transform_->inverse_transform_direction(h, work);
                H += (work);
            } else {
                if (!penalty_hessian(x, h)) {
                    return false;
                }

                apply_selection(h);

                if (has_scaling_matrix()) {
                    scaling_matrix()->apply(h, work);
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

        bool hessian(const Vector &x, Matrix &H) const {
            Vector work, h(layout(x), 0.);

            if (has_transform()) {
                // Transform to constraint base
                transform_->transform(x, work);
                if (!penalty_hessian(work, h)) {
                    return false;
                }

                apply_selection(h);

                // Scale contribuions
                if (has_scaling_matrix()) {
                    scaling_matrix()->apply(h, work);
                    h = work;
                }

                // Transform to problem base
                transform_->inverse_transform_direction(h, work);
                H.shift_diag(work);
            } else {
                if (!penalty_hessian(x, h)) {
                    return false;
                }

                apply_selection(h);

                if (has_scaling_matrix()) {
                    scaling_matrix()->apply(h, work);
                    H.shift_diag(work);
                } else {
                    H.shift_diag(h);
                }
            }

            return true;
        }

        /// FIXME shift is not considred in the energy value here
        bool value(const Vector &x, Scalar &value) const {
            Vector work, v(layout(x), 0.);

            if (has_transform()) {
                // Transform to constraint base
                transform_->transform(x, work);
                if (!penalty_value(work, v)) {
                    return false;
                }

                apply_selection(v);

                // Scale contribuions
                if (has_scaling_matrix()) {
                    scaling_matrix()->apply(v, work);
                    v = work;
                }

                // Transform to problem base
                transform_->inverse_transform_direction(v, work);
                value = sum(v);
            } else {
                if (!penalty_value(x, v)) {
                    return false;
                }

                apply_selection(v);

                if (has_scaling_matrix()) {
                    scaling_matrix()->apply(v, work);
                    value = sum(work);
                } else {
                    value = sum(v);
                }
            }

            return true;
        }

        bool gradient(const Vector &x, Vector &grad) const {
            Vector work, g(layout(x), 0.);

            if (has_transform()) {
                // Transform to constraint base
                transform_->transform(x, work);

                // Transform gradient
                transform_->transform(grad, g);

                if (!penalty_gradient(work, g)) {
                    return false;
                }

                apply_selection(g);

                // Scale contribuions
                if (has_scaling_matrix()) {
                    scaling_matrix()->apply(g, work);
                    g = work;
                }

                // Transform to problem base
                transform_->inverse_transform_direction(g, work);
                grad += work;
            } else {
                // Copy gradient
                g = grad;
                if (!penalty_gradient(x, g)) {
                    return false;
                }

                apply_selection(g);

                if (has_scaling_matrix()) {
                    scaling_matrix()->apply(g, work);
                    grad += work;
                } else {
                    grad += g;
                }
            }

            return true;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////

        inline bool has_transform() const { return static_cast<bool>(transform_); }
        void set_transform(const std::shared_ptr<Transformation> &t) { transform_ = t; }
        const std::shared_ptr<Transformation> &transform() const {
            assert(has_transform());
            return transform_;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////

        void set_selection(const std::shared_ptr<Vector> &selection) { selection_ = selection; }

        virtual void apply_selection(Vector &penalty_value) const {
            if (has_selection()) {
                penalty_value = e_mul(*selection_, penalty_value);
            }
        }

        bool has_selection() const { return static_cast<bool>(selection_); }

        inline const Vector &selection() const {
            assert(has_selection());
            return *selection_;
        }

        // void auto_selector(const bool val) { auto_selector_ = val; }

        void determine_boolean_selector() const {
            if (!this->box_) return;

            if (!selection_) {
                selection_ = std::make_shared<Vector>();
            }

            this->box_->determine_boolean_selector(-infinity_, infinity_, *selection_);
        }

        ////////////////////////////////////////////////////////////////////////////////////////////

        void set_box_constraints(const std::shared_ptr<BoxConstraints> &box) { box_ = box; }
        const std::shared_ptr<BoxConstraints> &box() { return box_; }

        void compute_diff_upper_bound(const Vector &x, Vector &diff) const {
            if (diff.empty()) {
                diff.zeros(layout(x));
            }

            auto ub_view = local_view_device(*box_->upper_bound());
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

            auto lb_view = local_view_device(*box_->lower_bound());
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

        inline bool verbose() const { return verbose_; }
        virtual void reset() {}

        inline bool has_scaling_matrix() const { return static_cast<bool>(scaling_matrix_); }

        inline void set_scaling_matrix(const std::shared_ptr<Matrix> &scaling_matrix) {
            scaling_matrix_ = scaling_matrix;
        }
        inline const std::shared_ptr<Matrix> &scaling_matrix() const {
            assert(has_scaling_matrix());
            return scaling_matrix_;
        }

        void read(Input &in) override {
            if (!Options()
                     .add_option(
                         "penalty_parameter", penalty_parameter_, "see Numerical Optimization - J. Nocedal, S. Wright.")
                     .add_option("verbose", verbose_, "Enable/Disable verbose output.")
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
            os << "verbose:\t" << verbose_ << "\n";

            os << "-----------------------------------------\n";
        }

        UTOPIA_NVCC_PRIVATE

        // Barrier requirements
        std::shared_ptr<BoxConstraints> box_;
        std::shared_ptr<Transformation> transform_;
        std::shared_ptr<Matrix> scaling_matrix_;
        std::shared_ptr<Vector> selection_;

        // Barrier parameters
        Scalar penalty_parameter_{1e-10};
        bool verbose_{false};
        Scalar infinity_{1e5};
    };
}  // namespace utopia

#endif  // UTOPIA_SHIFTED_PENALTY_HPP
