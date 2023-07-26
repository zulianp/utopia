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

        virtual bool hessian_and_gradient(const Vector &x, Matrix &H, Vector &g) const {
            return gradient(x, g) && hessian(x, H);
        }

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

        ////////////////////////////////////////////////////////////////////////////////////////////

        void set_box_constraints(const std::shared_ptr<BoxConstraints> &box) { box_ = box; }
        const std::shared_ptr<BoxConstraints> &box() const { return box_; }
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

        void read(Input &in) override {}

        void describe(std::ostream &os) const {}

        UTOPIA_NVCC_PRIVATE

        // Barrier requirements
        std::shared_ptr<BoxConstraints> box_;
        std::shared_ptr<Transformation> transform_;
        std::shared_ptr<Matrix> scaling_matrix_;
        std::shared_ptr<Vector> selection_;

        // Barrier parameters

        bool verbose_{false};
        Scalar infinity_{1e5};
    };
}  // namespace utopia

#endif  // UTOPIA_PENALTY_HPP
