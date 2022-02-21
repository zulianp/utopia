#ifndef UTOPIA_ACTIVE_SET_HPP
#define UTOPIA_ACTIVE_SET_HPP

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_BoxConstraints.hpp"
#include "utopia_DeviceView.hpp"

namespace utopia {
    template <class Vector>
    class ActiveSet : public Configurable {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using Layout = typename Traits<Vector>::Layout;

        virtual ~ActiveSet() = default;

        void expand_active_to_block(const int block_size) {
            auto d_indicator = local_view_device(indicator_);

            auto range = local_range_device(indicator_);
            RangeDevice<Vector> block_range(0, (range.end() - range.begin()) / block_size);

            assert(range.extent() == block_range.extent() * block_size);  // Check modulus

            parallel_for(
                block_range, UTOPIA_LAMBDA(const SizeType b) {
                    bool is_active = false;
                    for (int k = 0; k < block_size; ++k) {
                        SizeType i = b * block_size + k;
                        auto v = d_indicator.get(i);
                        if (v != 0) {
                            is_active = true;
                        }
                    }

                    if (is_active) {
                        for (int k = 0; k < block_size; ++k) {
                            SizeType i = b * block_size + k;
                            d_indicator.set(i, 1);
                        }
                    }
                });
        }

        virtual bool determine(const BoxConstraints<Vector> &c, const Vector &x) {
            if (!c.has_bound()) return false;

            if (!c.has_lower_bound()) {
                return determine_upper_bound(*c.upper_bound(), x);
            }

            auto d_lb = const_local_view_device(*c.lower_bound());
            auto d_ub = const_local_view_device(*c.upper_bound());
            auto d_x = const_local_view_device(x);

            auto d_indicator = local_view_device(indicator_);

            int changed = 0;
            parallel_reduce(
                local_range_device(indicator_),
                UTOPIA_LAMBDA(const SizeType i) {
                    const Scalar li = d_lb.get(i);
                    const Scalar ui = d_ub.get(i);
                    const Scalar xi = d_x.get(i);

                    const int old_ind = d_indicator.get(i);
                    const int new_ind = (xi <= (li + tol_) || xi >= (ui - tol_)) ? 1 : 0;

                    d_indicator.set(i, new_ind);

                    return old_ind != new_ind;
                },
                changed);

            changed = x.comm().sum(changed);

            if (verbose_ && changed && x.comm().rank() == 0) {
                utopia::out() << "ActiveSet changed: " << changed << std::endl;
            }

            return changed > 0;
        }

        void zero_out_active_rows(typename Traits<Vector>::Matrix &A, const Scalar value = 0.) const {
            set_zero_rows(A, indicator_, value);
        }

        void zero_out_active(Vector &x) const {
            auto d_i = const_local_view_device(indicator_);
            auto d_x = local_view_device(x);

            parallel_for(
                local_range_device(indicator_), UTOPIA_LAMBDA(const SizeType i) {
                    const Scalar is_active = d_i.get(i);
                    if (is_active == 1.0) {
                        d_x.set(i, 0.0);
                    }
                });
        }

        inline void init(const Layout &l) { indicator_.zeros(l); }
        inline const Vector &indicator() const { return indicator_; }
        inline Vector &indicator() { return indicator_; }
        inline void verbose(const bool verbose) { verbose_ = verbose; }
        inline void tol(const Scalar tol) { tol_ = tol; }
        inline SizeType num_active_sets() { return sum(indicator_); }

        void read(Input &in) override {
            in.get("verobse", verbose_);
            in.get("tol", tol_);
        }

        inline bool empty() const { return indicator_.empty(); }

    private:
        Vector indicator_;
        Scalar tol_{0};
        bool verbose_{false};

        bool determine_upper_bound(const Vector &upper_bound, const Vector &x) {
            auto d_ub = const_local_view_device(upper_bound);
            auto d_x = const_local_view_device(x);

            auto d_indicator = local_view_device(indicator_);

            int changed = 0;
            parallel_reduce(
                local_range_device(indicator_),
                UTOPIA_LAMBDA(const SizeType i) {
                    const Scalar ui = d_ub.get(i);
                    const Scalar xi = d_x.get(i);

                    const int old_ind = d_indicator.get(i);
                    const int new_ind = xi >= (ui - tol_) ? 1 : 0;

                    d_indicator.set(i, new_ind);

                    return old_ind != new_ind;
                },
                changed);

            changed = x.comm().sum(changed);

            if (verbose_ && changed && x.comm().rank() == 0) {
                utopia::out() << "ActiveSet changed: " << changed << std::endl;
            }

            return changed > 0;
        }
    };
}  // namespace utopia

#endif  // UTOPIA_ACTIVE_SET_HPP
