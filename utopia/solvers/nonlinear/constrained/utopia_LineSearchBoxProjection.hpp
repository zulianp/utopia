#ifndef UTOPIA_LINE_SEARCH_BOX_PROJECTION_HPP
#define UTOPIA_LINE_SEARCH_BOX_PROJECTION_HPP

#include "utopia_LS_Strategy.hpp"

#include <sstream>

namespace utopia {

    template <class Vector>
    class LineSearchBoxProjection : public LSStrategy<Vector> {
    public:
        using Super = utopia::LSStrategy<Vector>;
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        Scalar compute(const Vector &x, const Vector &correction) {
            if (!box_ || !box_->has_bound()) {
                assert(false);
                return 1;
            }

            Scalar alpha = 1;

            if (transform_) {
                Vector t_c;
                transform_->apply(correction, t_c);
                alpha = compute_aux(x, t_c);
            } else {
                alpha = compute_aux(x, correction);
            }

            if (this->verbose()) {
                std::stringstream ss;
                ss << "LineSearchBoxProjection::alpha = " << alpha << "\n";
                x.comm().root_print(ss.str());
            }

            return std::max(dumping_ * alpha, min_alpha_);
        }

        bool get_alpha(FunctionBase<Vector> &fun,
                       const Vector &g,
                       const Vector &x,
                       const Vector &correction,
                       Scalar &alpha) override {
            UTOPIA_UNUSED(fun);
            UTOPIA_UNUSED(g);

            alpha = compute(x, correction);
            return true;
        }

        bool get_alpha(LeastSquaresFunctionBase<Vector> &fun,
                       const Vector &g,
                       const Vector &x,
                       const Vector &correction,
                       Scalar &alpha) override {
            UTOPIA_UNUSED(fun);
            UTOPIA_UNUSED(g);

            alpha = compute(x, correction);
            return true;
        }

        void init_memory(const Layout &layout) override { buff_.zeros(layout); }

        void read(Input &in) override {
            Super::read(in);
            in.get("reach", reach_);
            in.get("dumping", dumping_);
            in.get("min_alpha", min_alpha_);
        }

        LineSearchBoxProjection(const std::shared_ptr<BoxConstraints<Vector>> &box,
                                const std::shared_ptr<Vector> &offset_vector = nullptr)
            : box_(box), offset_vector_(offset_vector) {
            // this->verbose(true);
        }

        inline void set_transform(const std::shared_ptr<Operator<Vector>> &transform) { transform_ = transform; }
        inline void set_box_constraints(const std::shared_ptr<BoxConstraints<Vector>> &box) { box_ = box; }
        inline void set_offset_vector(const std::shared_ptr<Vector> &offset_vector) { offset_vector_ = offset_vector; }

        /////////////
        UTOPIA_NVCC_PRIVATE
        /////////////

        // FIXME
        Vector buff_;
        Scalar reach_{1};
        Scalar dumping_{0.98};
        Scalar min_alpha_{1e-5};
        std::shared_ptr<BoxConstraints<Vector>> box_;

        // Usefull if we are solving for a delta_x
        std::shared_ptr<Vector> offset_vector_;
        std::shared_ptr<Operator<Vector>> transform_;

        void prepare_buff(const Vector &x) {
            if (transform_) {
                transform_->apply(x, buff_);
            } else {
                buff_ = x;
            }

            if (offset_vector_) {
                buff_ -= *offset_vector_;
            }
        }

        Scalar compute_aux(const Vector &x, const Vector &correction) {
            prepare_buff(x);

            if (box_->has_upper_bound() && box_->has_lower_bound()) {
                // Upper bound
                buff_ = (*box_->upper_bound()) - buff_;
                Scalar alpha_ub = compute_projection(x, correction);

                // Lower bound
                prepare_buff(x);
                buff_ -= (*box_->lower_bound());
                Scalar alpha_lb = compute_projection(x, correction);
                return std::min(alpha_lb, alpha_ub);

            } else if (box_->has_upper_bound()) {
                buff_ = (*box_->upper_bound()) - buff_;

                return compute_projection(x, correction);

            } else if (box_->has_lower_bound()) {
                buff_ -= (*box_->lower_bound());

                return compute_projection(x, correction);
            } else {
                return 1;
            }
        }

        Scalar compute_projection(const Vector &x, const Vector &correction) {
            UTOPIA_UNUSED(x);

            {
                auto correction_view = local_view_device(correction);
                auto buff_view = local_view_device(buff_);

                parallel_for(
                    local_range_device(correction), UTOPIA_LAMBDA(const SizeType i) {
                        auto xi = buff_view.get(i);
                        auto dxi = correction_view.get(i);

                        if (xi + dxi < (1 - reach_) * xi) {
                            auto val = -(xi * reach_) / dxi;
                            buff_view.set(i, val);
                        } else {
                            buff_view.set(i, 1);
                        }
                    });
            }

            return min(buff_);
        }
    };
}  // namespace utopia

#endif  // UTOPIA_LINE_SEARCH_BOX_PROJECTION_HPP
