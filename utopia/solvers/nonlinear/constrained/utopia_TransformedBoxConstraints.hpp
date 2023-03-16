#ifndef UTOPIA_TRANSFORMED_BOX_CONSTRAINTS_HPP
#define UTOPIA_TRANSFORMED_BOX_CONSTRAINTS_HPP

#include "utopia_BoxConstraints.hpp"
#include "utopia_Operator.hpp"

namespace utopia {
    template <class Matrix, class Vector>
    class Transformation {
    public:
        virtual ~Transformation() = default;
        virtual void transform(const Vector &in, Vector &out) = 0;
        virtual void transform(const Matrix &in, Matrix &out) = 0;

        virtual void transform_direction(const Vector &in, Vector &out) = 0;
        virtual void inverse_transform(const Vector &in, Vector &out) = 0;
        virtual void inverse_transform_direction(const Vector &in, Vector &out) = 0;
    };

    template <class Matrix, class Vector>
    class AffineTransformation : public Transformation<Matrix, Vector> {
    public:
        AffineTransformation(const std::shared_ptr<Matrix> &t,
                             const std::shared_ptr<Matrix> &t_inverse,
                             const std::shared_ptr<Vector> &offset_vector)
            : transform_(t), inverse_transform_(t_inverse), offset_vector_(offset_vector) {}

        void reset(const std::shared_ptr<Matrix> &t,
                   const std::shared_ptr<Matrix> &t_inverse,
                   const std::shared_ptr<Vector> &offset_vector) {
            transform_ = t;
            inverse_transform_ = t_inverse;
            offset_vector_ = offset_vector;
        }

        void transform_direction(const Vector &in, Vector &out) {
            if (transform_) {
                transform_->apply(in, out);
            } else {
                out = in;
            }
        }

        void inverse_transform_direction(const Vector &in, Vector &out) {
            if (transform_) {
                inverse_transform_->apply(in, out);
            } else {
                out = in;
            }
        }

        void apply_transform(const Matrix &in, Matrix &out) {
            if (transform_) {
                out = in;
            } else {
                out = transpose(*transform_) * out * *transform_;
            }
        }

        /// Transform to Box coordinates
        void apply_transform(const Vector &in, Vector &out) {
            if (!transform_ && !offset_vector_) {
                out = in;
                return;
            }

            if (offset_vector_) {
                work_ = in;
                work_ -= *offset_vector_;
                transform_direction(work_, out);
            } else {
                transform_direction(in, out);
            }
        }

        /// Transform to original coordinates
        void apply_inverse_transform(const Vector &in, Vector &out) {
            if (!inverse_transform_ && !offset_vector_) {
                out = in;
                return;
            }

            if (transform_) {
                inverse_transform_direction(in, out);
            } else {
                out = in;
            }

            if (offset_vector_) {
                out += *offset_vector_;
            }
        }

        inline const std::shared_ptr<Matrix> &transform() const { return transform_; }
        inline const std::shared_ptr<Vector> &offset_vector() const { return offset_vector_; }

    private:
        Vector work_;

        std::shared_ptr<Matrix> transform_;
        std::shared_ptr<Matrix> inverse_transform_;

        std::shared_ptr<Vector> offset_vector_;
    };

    template <class Matrix, class Vector>
    class TransformedBoxConstraints {
    public:
        void set_transform(const Transformation<Matrix, Vector> &t) { transform_ = t; }

        inline bool has_transform() const { return static_cast<bool>(transform_); }

        inline const std::shared_ptr<Transformation<Matrix, Vector>> &transform() const { return transform_; }

        inline void apply_selection(Vector &vec) const {
            if (selection_) {
                vec = e_mul(vec, *selection_);
            }
        }

        inline void set_selection(const std::shared_ptr<Vector> &selection) { selection_ = selection; }

        inline void set_box(const std::shared_ptr<BoxConstraints<Vector>> &box) { box_ = box; }
        inline const std::shared_ptr<BoxConstraints<Vector>> &box() const { return box_; }

    private:
        std::shared_ptr<BoxConstraints<Vector>> box_;
        std::shared_ptr<Transformation<Matrix, Vector>> transform_;
        std::shared_ptr<Vector> selection_;
    };
}  // namespace utopia

#endif  // UTOPIA_TRANSFORMED_BOX_CONSTRAINTS_HPP
