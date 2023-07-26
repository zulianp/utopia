#ifndef UTOPIA_LOG_BARRIER_FUNCTION_BASE_HPP
#define UTOPIA_LOG_BARRIER_FUNCTION_BASE_HPP

#include "utopia_BoxConstraints.hpp"
#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_Layout.hpp"
#include "utopia_Options.hpp"

#include "utopia_Penalty.hpp"
#include "utopia_TransformedBoxConstraints.hpp"

#include <limits>

namespace utopia {
    template <class Matrix, class Vector>
    class LogBarrierBase : public Penalty<Matrix, Vector> {
    public:
        using Super = utopia::Penalty<Matrix, Vector>;
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using BoxConstraints = utopia::BoxConstraints<Vector>;
        using Transformation = utopia::Transformation<Matrix, Vector>;

        LogBarrierBase() = default;
        explicit LogBarrierBase(const std::shared_ptr<BoxConstraints> &box) : Super(box) {}
        virtual ~LogBarrierBase() = default;

        virtual bool barrier_gradient(const Vector &x, Vector &g) const = 0;
        virtual bool barrier_hessian(const Vector &x, Vector &H) const = 0;
        virtual bool barrier_value(const Vector &x, Vector &v) const = 0;
        virtual std::string function_type() const = 0;

        virtual bool hessian_and_gradient(const Vector &x, Matrix &H, Vector &g) const final {
            return gradient(x, g) && hessian(x, H);
        }

        virtual bool hessian(const Vector &x, Vector &H) const final {
            Vector work, h(layout(x), 0.);

            if (this->has_transform()) {
                // Transform to constraint base
                this->transform()->transform(x, work);
                if (!barrier_hessian(work, h)) {
                    return false;
                }

                h *= current_barrier_parameter_;
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
                if (!barrier_hessian(x, h)) {
                    return false;
                }

                h *= current_barrier_parameter_;
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

        virtual bool hessian(const Vector &x, Matrix &H) const final {
            Vector work, h(layout(x), 0.);

            if (this->has_transform()) {
                // Transform to constraint base
                this->transform()->transform(x, work);
                if (!barrier_hessian(work, h)) {
                    return false;
                }

                h *= current_barrier_parameter_;

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
                if (!barrier_hessian(x, h)) {
                    return false;
                }

                h *= current_barrier_parameter_;

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

        virtual bool value(const Vector &x, Scalar &value) const final {
            Vector work, v(layout(x), 0.);

            if (this->has_transform()) {
                // Transform to constraint base
                this->transform()->transform(x, work);
                if (!barrier_value(work, v)) {
                    return false;
                }

                v *= current_barrier_parameter_;
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
                if (!barrier_value(x, v)) {
                    return false;
                }

                v *= current_barrier_parameter_;
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

        virtual bool gradient(const Vector &x, Vector &grad) const final {
            Vector work, g(layout(x), 0.);

            if (this->has_transform()) {
                // Transform to constraint base
                this->transform()->transform(x, work);
                if (!barrier_gradient(work, g)) {
                    return false;
                }

                g *= current_barrier_parameter_;
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
                if (!barrier_gradient(x, g)) {
                    return false;
                }

                g *= current_barrier_parameter_;
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

        virtual bool project_onto_feasibile_region(Vector &x) const {
            bool ok = false;
            if (this->has_transform()) {
                Vector temp_x;
                this->transform()->transform(x, temp_x);

                if (this->has_selection()) {
                    ok = extend_project_onto_feasibile_region_with_selection(temp_x);
                } else {
                    ok = extend_project_onto_feasibile_region(temp_x);
                }

                this->transform()->inverse_transform(temp_x, x);

            } else {
                if (this->has_selection()) {
                    ok = extend_project_onto_feasibile_region_with_selection(x);
                } else {
                    ok = extend_project_onto_feasibile_region(x);
                }
            }

            return ok;
        }

        bool extend_project_onto_feasibile_region_with_selection(Vector &x) const {
            if (this->box()->has_upper_bound()) {
                auto ub_view = local_view_device(*this->box()->upper_bound());
                auto x_view = local_view_device(x);
                auto selector_view = local_view_device(this->selection());

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

            if (this->box()->has_lower_bound()) {
                auto lb_view = local_view_device(*this->box()->lower_bound());
                auto x_view = local_view_device(x);
                auto selector_view = local_view_device(this->selection());

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

        bool extend_project_onto_feasibile_region(Vector &x) const {
            if (this->box()->has_upper_bound()) {
                auto ub_view = local_view_device(*this->box()->upper_bound());
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

            if (this->box()->has_lower_bound()) {
                auto lb_view = local_view_device(*this->box()->lower_bound());
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

        ////////////////////////////////////////////////////////////////////////////////////////////

        // inline bool this->has_transform() const { return static_cast<bool>(transform_); }
        // void set_transform(const std::shared_ptr<Transformation> &t) { transform_ = t; }
        // const std::shared_ptr<Transformation> &transform() const {
        //     assert(this->has_transform());
        //     return transform_;
        // }

        ////////////////////////////////////////////////////////////////////////////////////////////

        // void set_selection(const std::shared_ptr<Vector> &selection) { selection_ = selection; }

        // virtual void this->apply_selection(Vector &barrier_value) const {
        //     if (this->has_selection()) {
        //         barrier_value = e_mul(*selection_, barrier_value);
        //     }
        // }

        // bool this->has_selection() const { return static_cast<bool>(selection_); }

        // inline const Vector &selection() const {
        //     assert(this->has_selection());
        //     return *selection_;
        // }

        // void auto_selector(const bool val) { auto_selector_ = val; }

        // void determine_boolean_selector() const {
        //     if (!this->box_) return;

        //     if (!selection_) {
        //         selection_ = std::make_shared<Vector>();
        //     }

        //     this->box()->determine_boolean_selector(-infinity_, infinity_, *selection_);
        // }

        ////////////////////////////////////////////////////////////////////////////////////////////

        // void set_box_constraints(const std::shared_ptr<BoxConstraints> &box) { box_ = box; }
        // const std::shared_ptr<BoxConstraints> &box() { return box_; }

        virtual void update_barrier() {
            current_barrier_parameter_ =
                std::max(current_barrier_parameter_ * barrier_parameter_shrinking_factor_, min_barrier_parameter_);
        }

        void compute_diff_upper_bound(const Vector &x, Vector &diff) const {
            if (diff.empty()) {
                diff.zeros(layout(x));
            }

            auto ub_view = local_view_device(*this->box()->upper_bound());
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

                    if (d < 0) {
                        if (this->verbose()) printf("Negative distance at dof %d: %g, ILLEGAL!\n", int(i), d);
                    }

                    // assert(d > 0);

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

            auto zero = zero_;
            parallel_for(
                local_range_device(x), UTOPIA_LAMBDA(const SizeType i) {
                    auto xi = x_view.get(i);
                    auto lbi = lb_view.get(i);
                    auto d = xi - lbi;
                    if (d == 0.) {
                        d = zero;
                    }

                    // assert(d > 0);

                    diff_view.set(i, d);
                });
        }

        // inline bool verbose() const { return verbose_; }
        // virtual void reset() { current_barrier_parameter_ = barrier_parameter_; }

        // inline bool this->has_scaling_matrix() const { return static_cast<bool>(scaling_matrix_); }

        // inline void set_scaling_matrix(const std::shared_ptr<Matrix> &scaling_matrix) {
        //     scaling_matrix_ = scaling_matrix;
        // }
        // inline const std::shared_ptr<Matrix> &scaling_matrix() const {
        //     assert(this->has_scaling_matrix());
        //     return scaling_matrix_;
        // }

        void read(Input &in) override {
            Super::read(in);

            if (!Options()
                     .add_option(
                         "barrier_parameter", barrier_parameter_, "see Numerical Optimization - J. Nocedal, S. Wright.")
                     .add_option("barrier_parameter_shrinking_factor",
                                 barrier_parameter_shrinking_factor_,
                                 "Factor with which the barrier term is reduced.")
                     .add_option("min_barrier_parameter", min_barrier_parameter_, "Smallest barrier parameter allowed.")
                     .add_option(
                         "soft_boundary", soft_boundary_, "A value in (0,0.01) used to project in the feasible region.")
                     // .add_option("verbose", verbose_, "Enable/Disable verbose output.")
                     .add_option("zero", zero_, "numerical zero (>0)")
                     // .add_option(
                     // "auto_selector", auto_selector_, "Determine selection automatically (use it for testing only)")
                     .parse(in)) {
                return;
            }

            if (mpi_world_rank() == 0) {
                describe(utopia::out().stream());
            }
        }

        void describe(std::ostream &os) const {
            os << "-----------------------------------------\n";
            os << "utopia::LogBarrierBase\n";
            os << "-----------------------------------------\n";

            os << "barrier_parameter:\t" << barrier_parameter_ << "\n";
            os << "barrier_parameter_shrinking_factor:\t" << barrier_parameter_shrinking_factor_ << "\n";
            os << "min_barrier_parameter:\t" << min_barrier_parameter_ << "\n";
            os << "current_barrier_parameter:\t" << current_barrier_parameter_ << "\n";

            // Soft limit
            os << "soft_boundary:\t" << soft_boundary_ << "\n";
            os << "zero:\t" << zero_ << "\n";

            // Auto selector tol
            // os << "infinity:\t" << infinity_ << "\n";
            // os << "auto_selector:\t" << auto_selector_ << "\n";

            // os << "verbose:\t" << verbose_ << "\n";

            os << "-----------------------------------------\n";
        }

        void set_all_barrier_parameters_to(const Scalar value) {
            barrier_parameter_ = value;
            min_barrier_parameter_ = value;
            current_barrier_parameter_ = value;
        }

        UTOPIA_NVCC_PRIVATE

        // Barrier requirements
        // std::shared_ptr<LogBarrierBase> barrier_;
        // std::shared_ptr<BoxConstraints> box_;
        // std::shared_ptr<Transformation> transform_;
        // std::shared_ptr<Matrix> scaling_matrix_;
        // std::shared_ptr<Vector> selection_;

        // Barrier parameters
        Scalar barrier_parameter_{1e-10};
        Scalar barrier_parameter_shrinking_factor_{0.1};
        Scalar min_barrier_parameter_{1e-10};
        Scalar current_barrier_parameter_{1e-10};

        // Soft limit
        Scalar soft_boundary_{1e-7};
        Scalar zero_{1e-20};

        // // Auto selector tol
        // Scalar infinity_{std::numeric_limits<Scalar>::max()};
        // bool auto_selector_{false};

        // bool verbose_{false};
    };

}  // namespace utopia

#endif  // UTOPIA_LOG_BARRIER_FUNCTION_BASE_HPP
