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
    class BoundedLogBarrier : public LogBarrierBase<Matrix, Vector> {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Function = utopia::Function<Matrix, Vector>;
        using BoxConstraints = utopia::BoxConstraints<Vector>;
        using Super = utopia::LogBarrierBase<Matrix, Vector>;

        class BarrierIPC {
        public:
            UTOPIA_INLINE_FUNCTION Scalar value(const Scalar d, Scalar mu) const {
                return -(d < d_hat) * (mu * (d - d_hat) * (d - d_hat) * device::log(d / d_hat));
            }

            UTOPIA_INLINE_FUNCTION Scalar gradient(const Scalar d, Scalar mu) const {
                return -(d < d_hat) * (mu * (d - d_hat) * ((d - d_hat) + 2 * d * device::log(d / d_hat)) / d);
            }
            UTOPIA_INLINE_FUNCTION Scalar hessian(const Scalar d, Scalar mu) const {
                return -(d < d_hat) * mu *
                       (-(d_hat * d_hat) / (d * d) - 2 * d_hat / d + 2 * device::log(d / d_hat) + 3);
            }

            Scalar d_hat{0.1};
        };

        class BarrierMine {
        public:
            UTOPIA_INLINE_FUNCTION Scalar value(const Scalar d, Scalar mu) const {
                auto temp = ((d - d_hat) / d_hat);
                return -(d < d_hat) * (mu * temp * temp * log(d / d_hat));
            }

            UTOPIA_INLINE_FUNCTION Scalar gradient(const Scalar d, Scalar mu) const {
                return -(d < d_hat) *
                       (mu * (d - d_hat) * ((d - d_hat) + 2 * d * device::log(d / d_hat)) / (d * d_hat * d_hat));
            }

            UTOPIA_INLINE_FUNCTION Scalar hessian(const Scalar d, Scalar mu) const {
                auto d2 = d * d;
                auto d_hat2 = d_hat * d_hat;
                return -(d < d_hat) * mu * (2 * d2 * device::log(d / d_hat) + 3 * d2 - 2 * d * d_hat - d_hat2) /
                       (d2 * d_hat2);
            }

            Scalar d_hat{0.1};
        };

        // using DefaultBarrier = BarrierIPC;
        using DefaultBarrier = BarrierMine;

        BoundedLogBarrier() = default;
        explicit BoundedLogBarrier(const std::shared_ptr<BoxConstraints> &box) : Super(box) {}

        void reset() override {}

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

        void add_barrier_value(const Vector &diff, Scalar &val) const {
            auto diff_view = local_view_device(diff);

            auto d_hat = barrier_thickness_;

            // Currently it is not adaptive like in the paper
            auto stiffness = barrier_stiffness();

            DefaultBarrier b{d_hat};

            Scalar b_val = 0.;
            parallel_reduce(
                local_range_device(diff),
                UTOPIA_LAMBDA(const SizeType i) {
                    auto d_i = diff_view.get(i);
                    return b.value(d_i, stiffness);
                },
                b_val);

            val += b_val;
        }

        void add_barrier_gradient(const Vector &diff, Vector &g) const {
            auto diff_view = local_view_device(diff);
            auto g_view = local_view_device(g);

            auto d_hat = barrier_thickness_;

            // Currently it is not adaptive like in the paper
            auto stiffness = barrier_stiffness();

            DefaultBarrier b{d_hat};

            parallel_for(
                local_range_device(diff), UTOPIA_LAMBDA(const SizeType i) {
                    auto d_i = diff_view.get(i);

                    auto b_g = b.gradient(d_i, stiffness);

                    assert(b_g == b_g);

                    auto g_i = g_view.get(i);
                    g_view.set(i, g_i - b_g);
                });
        }

        // !!! diff is modified inside !!!
        void add_barrier_hessian(Vector &diff, Matrix &hessian) const {
            {
                auto diff_view = local_view_device(diff);
                auto d_hat = barrier_thickness_;

                // Currently it is not adaptive like in the paper
                auto stiffness = barrier_stiffness();

                DefaultBarrier b{d_hat};

                parallel_for(
                    local_range_device(diff), UTOPIA_LAMBDA(const SizeType i) {
                        auto d_i = diff_view.get(i);
                        auto b_H = b.hessian(d_i, stiffness);

                        assert(b_H == b_H);

                        diff_view.set(i, b_H);
                    });
            }

            hessian.shift_diag(diff);
        }

        void hessian_and_gradient(const Vector &x, Matrix &H, Vector &g) const override {
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

        void hessian(const Vector &x, Matrix &H) const override {
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

        void gradient(const Vector &x, Vector &g) const override {
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

        void value(const Vector &x, Scalar &value) const override {
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

        inline Scalar barrier_stiffness() const {
            return this->current_barrier_parameter_;
            // return 1;
        }

    private:
        Scalar barrier_thickness_{0.01};
    };

    template <class Matrix, class Vector>
    class BoundedLogBarrierFunction : public LogBarrierFunctionBase<Matrix, Vector> {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Function = utopia::Function<Matrix, Vector>;
        using BoxConstraints = utopia::BoxConstraints<Vector>;
        using Super = utopia::LogBarrierFunctionBase<Matrix, Vector>;
        using BoundedLogBarrier = utopia::BoundedLogBarrier<Matrix, Vector>;

        BoundedLogBarrierFunction() { this->set_barrier(std::make_shared<BoundedLogBarrier>()); }

        BoundedLogBarrierFunction(const std::shared_ptr<Function> &unconstrained,
                                  const std::shared_ptr<BoxConstraints> &box)
            : Super(unconstrained, std::make_shared<BoundedLogBarrier>(box)) {}

        inline std::string function_type() const override { return "BoundedLogBarrierFunction"; }
    };

}  // namespace utopia

#endif  // UTOPIA_BOUNDED_LOG_BARRIER_FUNCTION_HPP
