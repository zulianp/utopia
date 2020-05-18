#ifndef UTOPIA_BOX_CONSTRAINTS_HPP
#define UTOPIA_BOX_CONSTRAINTS_HPP

#include "utopia_Base.hpp"

#include "utopia_Clonable.hpp"
#include "utopia_Factory.hpp"
#include "utopia_Layout.hpp"
#include "utopia_Traits.hpp"

#include <limits>
#include <memory>
#include <utility>

namespace utopia {

    template <class Vector>
    class BoxConstraints : public Clonable {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        BoxConstraints(std::shared_ptr<Vector> lower_bound, std::shared_ptr<Vector> upper_bound)
            : lower_bound_(std::move(lower_bound)),
              upper_bound_(std::move(upper_bound)),
              min_val_(-std::numeric_limits<Scalar>::max()),
              max_val_(std::numeric_limits<Scalar>::max()) {}

        BoxConstraints()
            : min_val_(-std::numeric_limits<Scalar>::max()), max_val_(std::numeric_limits<Scalar>::max()) {}

        BoxConstraints(const BoxConstraints &other)
            : min_val_(other.min_val_), max_val_(other.max_val_), uniform_(other.uniform_) {
            if (other.lower_bound_) {
                lower_bound_ = std::make_shared<Vector>(*other.lower_bound_);
            }

            if (other.upper_bound_) {
                upper_bound_ = std::make_shared<Vector>(*other.upper_bound_);
            }
        }

        void set(const std::shared_ptr<Vector> &lower_bound, const std::shared_ptr<Vector> &upper_bound) {
            lower_bound_ = lower_bound;
            upper_bound_ = upper_bound;
        }

        BoxConstraints(BoxConstraints &&other)
            : min_val_(other.min_val_), max_val_(other.max_val_), uniform_(other.uniform_) {
            if (other.lower_bound_) {
                lower_bound_ = std::move(other.lower_bound_);
            }

            if (other.upper_bound_) {
                upper_bound_ = std::move(other.upper_bound_);
            }
        }

        BoxConstraints &operator=(const BoxConstraints &other) {
            if (this == &other) {
                return *this;
            }

            min_val_ = other.min_val_;
            max_val_ = other.max_val_;
            uniform_ = other.uniform_;
            if (other.lower_bound_) {
                lower_bound_ = std::make_shared<Vector>(*other.lower_bound_);
            }

            if (other.upper_bound_) {
                upper_bound_ = std::make_shared<Vector>(*other.upper_bound_);
            }

            return *this;
        }

        BoxConstraints &operator=(BoxConstraints &&other) {
            min_val_ = other.min_val_;
            max_val_ = other.max_val_;
            uniform_ = other.uniform_;
            if (other.lower_bound_) {
                lower_bound_ = std::move(other.lower_bound_);
            }

            if (other.upper_bound_) {
                upper_bound_ = std::move(other.upper_bound_);
            }

            return *this;
        }

        BoxConstraints *clone() const override { return new BoxConstraints(*this); }

        ~BoxConstraints() override = default;

        inline std::shared_ptr<Vector> &upper_bound() { return upper_bound_; }

        inline std::shared_ptr<const Vector> upper_bound() const { return upper_bound_; }

        inline std::shared_ptr<Vector> &lower_bound() { return lower_bound_; }

        inline std::shared_ptr<const Vector> lower_bound() const { return lower_bound_; }

        inline bool has_lower_bound() const { return static_cast<bool>(lower_bound_); }

        inline bool has_upper_bound() const { return static_cast<bool>(upper_bound_); }

        inline bool has_bound() const { return has_lower_bound() || has_upper_bound(); }

        inline bool has_bounds() const { return has_lower_bound() && has_upper_bound(); }

        inline bool has_empty_bounds() const { return !has_lower_bound() || !has_upper_bound(); }

        inline bool valid(const Layout &lo) {
            return !has_empty_bounds() && lo.same(layout(*lower_bound())) && lo.same(layout(*upper_bound()));
        }

        inline void fill_empty_bounds(const Layout &lo) {
            if (this->has_bounds()) {
                if (!lo.same(layout(*lower_bound_))) {
                    lower_bound_ = std::make_shared<Vector>();
                    lower_bound_->values(lo, min_val_);
                }

                // couldn't this 2nd check be avoided?
                if (!lo.same(layout(*upper_bound_))) {
                    upper_bound_ = std::make_shared<Vector>();
                    upper_bound_->values(lo, max_val_);
                }

            } else if (!lower_bound_ && !upper_bound_) {
                lower_bound_ = std::make_shared<Vector>();
                upper_bound_ = std::make_shared<Vector>();

                lower_bound_->values(lo, min_val_);
                upper_bound_->values(lo, max_val_);

            } else {
                if (!lower_bound_) {
                    lower_bound_ = std::make_shared<Vector>();
                    lower_bound_->values(lo, min_val_);
                }

                if (!upper_bound_) {
                    upper_bound_ = std::make_shared<Vector>();
                    upper_bound_->values(lo, max_val_);
                }
            }
        }

        inline bool uniform() const { return uniform_; }

        inline void uniform(const bool &flg) { uniform_ = flg; }

    private:
        std::shared_ptr<Vector> lower_bound_;
        std::shared_ptr<Vector> upper_bound_;
        Scalar min_val_ = -9e9;
        Scalar max_val_ = 9e9;
        bool uniform_{false};
    };

    template <class Vector>
    inline BoxConstraints<Vector> make_box_constaints(const std::shared_ptr<Vector> &lower_bound,
                                                      const std::shared_ptr<Vector> &upper_bound) {
        return BoxConstraints<Vector>(lower_bound, upper_bound);
    }

    template <class Vector>
    inline BoxConstraints<Vector> make_lower_bound_constraints(const std::shared_ptr<Vector> &lower_bound) {
        return BoxConstraints<Vector>(lower_bound, nullptr);
    }

    template <class Vector>
    inline BoxConstraints<Vector> make_upper_bound_constraints(const std::shared_ptr<Vector> &upper_bound) {
        return BoxConstraints<Vector>(nullptr, upper_bound);
    }

}  // namespace utopia

#endif  // UTOPIA_BOX_CONSTRAINTS_HPP
