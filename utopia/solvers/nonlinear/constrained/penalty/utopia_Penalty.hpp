#ifndef UTOPIA_PENALTY_HPP
#define UTOPIA_PENALTY_HPP

#include "utopia_BoxConstraints.hpp"
#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_Layout.hpp"
#include "utopia_Options.hpp"
#include "utopia_TransformedBoxConstraints.hpp"

#include <limits>

namespace utopia {

    template <class Matrix, class Vector = typename Traits<Matrix>::Vector>
    class Penalty : public Configurable {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using BoxConstraints = utopia::BoxConstraints<Vector>;
        using Transformation = utopia::Transformation<Matrix, Vector>;

        Penalty() = default;
        explicit Penalty(const std::shared_ptr<BoxConstraints> &box) : box_(box) {}
        virtual ~Penalty() = default;

        virtual std::string function_type() const = 0;

        virtual bool hessian_and_gradient(const Vector &x, Matrix &H, Vector &g) const {
            return gradient(x, g) && hessian(x, H);
        }

        virtual bool hessian(const Vector &x, Vector &diagH) const = 0;
        virtual bool hessian(const Vector &x, Matrix &H) const = 0;
        virtual bool value(const Vector &x, Scalar &value) const = 0;
        virtual bool gradient(const Vector &x, Vector &grad) const = 0;

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

        ////////////////////////////////////////////////////////////////////////////////////////////
        virtual void set_auxiliary_forcing(const std::shared_ptr<Vector> & /*vec*/) {}
        virtual bool supports_auxiliary_forcing() const { return false; }

        void set_box_constraints(const std::shared_ptr<BoxConstraints> &box) { box_ = box; }
        const std::shared_ptr<BoxConstraints> &box() const { return box_; }
        inline bool verbose() const { return verbose_; }
        virtual void reset() {}

        virtual void update() {}

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
                         "soft_boundary", soft_boundary_, "A value in (0,0.01) used to project in the feasible region.")
                     .parse(in)) {
                return;
            }
        }

        void describe(std::ostream &os) const { os << "soft_boundary:\t" << soft_boundary_ << "\n"; }

        UTOPIA_NVCC_PRIVATE

        // Barrier requirements
        std::shared_ptr<BoxConstraints> box_;
        std::shared_ptr<Transformation> transform_;
        std::shared_ptr<Matrix> scaling_matrix_;
        std::shared_ptr<Vector> selection_;

        // Barrier parameters

        bool verbose_{false};
        Scalar infinity_{1e5};
        Scalar soft_boundary_{1e-7};
    };
}  // namespace utopia

#endif  // UTOPIA_PENALTY_HPP
