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
            UTOPIA_INLINE_FUNCTION Scalar value(const Scalar d) const {
                return -(d < d_hat) * ((d - d_hat) * (d - d_hat) * device::log(d / d_hat));
            }

            UTOPIA_INLINE_FUNCTION Scalar gradient(const Scalar d) const {
                return -(d < d_hat) * ((d - d_hat) * ((d - d_hat) + 2 * d * device::log(d / d_hat)) / d);
            }

            UTOPIA_INLINE_FUNCTION Scalar hessian(const Scalar d) const {
                return -(d < d_hat) * (-(d_hat * d_hat) / (d * d) - 2 * d_hat / d + 2 * device::log(d / d_hat) + 3);
            }

            Scalar d_hat{0.1};
        };

        class BarrierMine {
        public:
            UTOPIA_INLINE_FUNCTION Scalar value(const Scalar d) const {
                auto temp = ((d - d_hat) / d_hat);
                return -(d < d_hat) * (temp * temp * log(d / d_hat));
            }

            UTOPIA_INLINE_FUNCTION Scalar gradient(const Scalar d) const {
                return -(d < d_hat) *
                       ((d - d_hat) * ((d - d_hat) + 2 * d * device::log(d / d_hat)) / (d * d_hat * d_hat));
            }

            UTOPIA_INLINE_FUNCTION Scalar hessian(const Scalar d) const {
                auto d2 = d * d;
                auto d_hat2 = d_hat * d_hat;
                return -(d < d_hat) * (2 * d2 * device::log(d / d_hat) + 3 * d2 - 2 * d * d_hat - d_hat2) /
                       (d2 * d_hat2);
            }

            Scalar d_hat{0.1};
        };

        class PolynomialBarrier {
        public:
            UTOPIA_INLINE_FUNCTION Scalar value(const Scalar d) const {
                auto param = (d - d_hat) / d_hat;
                return (d < d_hat) * 0.5 * (param * param);
            }

            UTOPIA_INLINE_FUNCTION Scalar gradient(const Scalar d) const {
                auto param = (d - d_hat) / d_hat;
                return (d < d_hat) * (param);
            }

            UTOPIA_INLINE_FUNCTION Scalar hessian(const Scalar d) const { return (d < d_hat); }

            Scalar d_hat{0.1};
        };

        // using DefaultBarrier = BarrierIPC;
        using DefaultBarrier = BarrierMine;

        BoundedLogBarrier() = default;
        explicit BoundedLogBarrier(const std::shared_ptr<BoxConstraints> &box) : Super(box) {}

        void reset() override { Super::reset(); }

        void read(Input &in) override {
            Super::read(in);
            if (!Options()
                     .add_option("barrier_thickness",
                                 barrier_thickness_,
                                 "see: Technical Supplement to Incremental Potential Contact: Intersection- and "
                                 "Inversion-free, Large-Deformation Dynamics.")
                     .add_option("barrier_subtype", barrier_subtype, "{default|polynomial}")
                     .parse(in)) {
                return;
            }
        }

        inline std::string function_type() const override { return "BoundedLogBarrier"; }

        bool barrier_hessian(const Vector &x, Vector &h) const override {
            if (barrier_subtype == "polynomial") {
                PolynomialBarrier b{barrier_thickness_};
                return barrier_hessian_aux(b, x, h);
            } else {
                DefaultBarrier b{barrier_thickness_};
                return barrier_hessian_aux(b, x, h);
            }
        }

        bool barrier_gradient(const Vector &x, Vector &g) const override {
            if (barrier_subtype == "polynomial") {
                PolynomialBarrier b{barrier_thickness_};
                return barrier_gradient_aux(b, x, g);
            } else {
                DefaultBarrier b{barrier_thickness_};
                return barrier_gradient_aux(b, x, g);
            }
        }

        bool barrier_value(const Vector &x, Vector &value) const override {
            if (barrier_subtype == "polynomial") {
                PolynomialBarrier b{barrier_thickness_};
                return barrier_value_aux(b, x, value);
            } else {
                DefaultBarrier b{barrier_thickness_};
                return barrier_value_aux(b, x, value);
            }
        }

        UTOPIA_NVCC_PRIVATE
        Scalar barrier_thickness_{0.01};
        std::string barrier_subtype{"default"};

        //////////////////////// Helper methods ////////////////////////

        template <class Barrier>
        void in_place_barrier_value(Vector &diff_in_value_out, Barrier b) const {
            auto view = local_view_device(diff_in_value_out);

            parallel_for(
                local_range_device(diff_in_value_out), UTOPIA_LAMBDA(const SizeType i) {
                    auto d_i = view.get(i);
                    auto b_v = b.value(d_i);

                    assert(b_v == b_v);  // NaN check

                    view.set(i, b_v);
                });
        }

        template <class Barrier>
        void in_place_barrier_gradient(Vector &diff_in_gradient_out, Barrier b) const {
            auto view = local_view_device(diff_in_gradient_out);

            parallel_for(
                local_range_device(diff_in_gradient_out), UTOPIA_LAMBDA(const SizeType i) {
                    auto d_i = view.get(i);
                    auto b_g = b.gradient(d_i);

                    assert(b_g == b_g);  // NaN check

                    view.set(i, -b_g);
                });
        }

        template <class Barrier>
        void in_place_barrier_hessian(Vector &diff_in_hessian_out, Barrier b) const {
            auto view = local_view_device(diff_in_hessian_out);

            parallel_for(
                local_range_device(diff_in_hessian_out), UTOPIA_LAMBDA(const SizeType i) {
                    auto d_i = view.get(i);
                    auto b_H = b.hessian(d_i);

                    assert(b_H == b_H);  // NaN check

                    view.set(i, b_H);
                });
        }

        /////////////////////////////////////////

        template <class Barrier>
        bool barrier_hessian_aux(Barrier b, const Vector &x, Vector &h) const {
            Vector work;

            if (this->box_->has_upper_bound()) {
                this->compute_diff_upper_bound(x, work);

                in_place_barrier_hessian(work, b);

                if (h.empty()) {
                    h.zeros(layout(work));
                }

                h += work;
            }

            if (this->box_->has_lower_bound()) {
                this->compute_diff_lower_bound(x, work);

                in_place_barrier_hessian(work, b);

                if (h.empty()) {
                    h.zeros(layout(work));
                }

                h += work;
            }

            return true;
        }

        template <class Barrier>
        bool barrier_gradient_aux(Barrier b, const Vector &x, Vector &g) const {
            Vector work;

            if (this->box_->has_upper_bound()) {
                this->compute_diff_upper_bound(x, work);

                in_place_barrier_gradient(work, b);

                g += work;
            }

            if (this->box_->has_lower_bound()) {
                this->compute_diff_lower_bound(x, work);

                in_place_barrier_gradient(work, b);

                g -= work;
            }

            return true;
        }

        template <class Barrier>
        bool barrier_value_aux(Barrier b, const Vector &x, Vector &value) const {
            Vector work;

            if (this->box_->has_upper_bound()) {
                work = *this->box_->upper_bound() - x;

                in_place_barrier_value(work, b);

                value += work;
            }

            if (this->box_->has_lower_bound()) {
                work = x - *this->box_->lower_bound();

                in_place_barrier_value(work, b);

                value += work;
            }

            return true;
        }
    };

}  // namespace utopia

#endif  // UTOPIA_BOUNDED_LOG_BARRIER_FUNCTION_HPP
