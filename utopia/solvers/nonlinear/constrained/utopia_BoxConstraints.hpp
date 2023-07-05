#ifndef UTOPIA_BOX_CONSTRAINTS_HPP
#define UTOPIA_BOX_CONSTRAINTS_HPP

#include "utopia_Base.hpp"

#include "utopia_Clonable.hpp"
#include "utopia_Factory.hpp"
#include "utopia_Layout.hpp"
#include "utopia_Traits.hpp"

#include <limits>
#include <memory>
#include <sstream>
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

        BoxConstraints(const Scalar &min_value, const Scalar &max_value, const Layout &lo)
            : min_val_(min_value), max_val_(max_value), uniform_(true) {
            lower_bound_ = std::make_shared<Vector>();
            lower_bound_->values(lo, min_value);

            upper_bound_ = std::make_shared<Vector>();
            upper_bound_->values(lo, max_value);
        }

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

        SizeType local_count_active(const Vector &x, const Scalar tol = 0) const {
            SizeType count_active = 0;

            auto x_view = const_local_view_device(x);
            auto r = local_range_device(x);

            if (this->has_upper_bound() && this->has_lower_bound()) {
                auto u_view = const_local_view_device(*this->upper_bound());
                auto l_view = const_local_view_device(*this->lower_bound());

                parallel_reduce(
                    r,
                    UTOPIA_LAMBDA(const SizeType &i)->SizeType {
                        return (l_view.get(i) + tol >= x_view.get(i)) || (u_view.get(i) - tol <= x_view.get(i));
                    },
                    count_active);

            } else if (this->has_upper_bound()) {
                auto u_view = const_local_view_device(*this->upper_bound());

                parallel_reduce(
                    r,
                    UTOPIA_LAMBDA(const SizeType &i)->SizeType { return (u_view.get(i) - tol <= x_view.get(i)); },
                    count_active);

            } else if (this->has_lower_bound()) {
                auto l_view = const_local_view_device(*this->lower_bound());

                parallel_reduce(
                    r,
                    UTOPIA_LAMBDA(const SizeType &i)->SizeType { return (l_view.get(i) + tol >= x_view.get(i)); },
                    count_active);
            }

            return count_active;
        }

        SizeType count_active(const Vector &x, const Scalar tol = 0) const {
            SizeType ca = local_count_active(x, tol);
            return x.comm().sum(ca);
        }

        bool determine_boolean_selector(const Scalar &negative_infinity,
                                        const Scalar &positive_infinity,
                                        Vector &boolean_selector,
                                        const bool verbose = false) const {
            if (!has_bound()) return false;

            if (empty(boolean_selector)) {
                if (this->has_upper_bound()) {
                    boolean_selector.zeros(layout(*this->upper_bound()));
                } else if (this->has_lower_bound()) {
                    boolean_selector.zeros(layout(*this->lower_bound()));
                }
            }

            if (this->has_upper_bound() && this->has_lower_bound()) {
                auto lb_view = local_view_device(*this->lower_bound());
                auto ub_view = local_view_device(*this->upper_bound());
                auto selector_view = local_view_device(boolean_selector);

                parallel_for(
                    local_range_device(boolean_selector), UTOPIA_LAMBDA(const SizeType i) {
                        auto ubi = ub_view.get(i);
                        auto lbi = lb_view.get(i);

                        if (ubi < positive_infinity || lbi > negative_infinity) {
                            selector_view.set(i, 1.);
                        }
                    });

            } else if (this->has_upper_bound()) {
                auto ub_view = local_view_device(*this->upper_bound());
                auto selector_view = local_view_device(boolean_selector);

                parallel_for(
                    local_range_device(boolean_selector), UTOPIA_LAMBDA(const SizeType i) {
                        auto bound = ub_view.get(i);

                        if (bound < positive_infinity) {
                            selector_view.set(i, 1.);
                        }
                    });
            } else if (this->has_lower_bound()) {
                auto lb_view = local_view_device(*this->lower_bound());
                auto selector_view = local_view_device(boolean_selector);

                parallel_for(
                    local_range_device(boolean_selector), UTOPIA_LAMBDA(const SizeType i) {
                        auto bound = lb_view.get(i);

                        if (bound > negative_infinity) {
                            selector_view.set(i, 1.);
                        }
                    });
            }

            if (verbose) {
                Scalar sum_selected = sum(boolean_selector);
                if (boolean_selector.comm().rank() == 0) {
                    std::stringstream ss;
                    ss << "Selected: " << SizeType(sum_selected) << "\n";
                    boolean_selector.comm().root_print(ss.str());
                }
            }

            return true;
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

        inline Scalar min_value() { return min_val_; }
        inline Scalar max_value() { return max_val_; }

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

    template <class Vector>
    inline BoxConstraints<Vector> make_uniform_box_constraints(const typename Traits<Vector>::Scalar &min_value,
                                                               const typename Traits<Vector>::Scalar &max_value,
                                                               const typename Traits<Vector>::Layout &lo) {
        return BoxConstraints<Vector>(min_value, max_value, lo);
    }

}  // namespace utopia

#endif  // UTOPIA_BOX_CONSTRAINTS_HPP
