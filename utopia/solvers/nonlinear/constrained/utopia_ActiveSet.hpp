#ifndef UTOPIA_ACTIVE_SET_HPP
#define UTOPIA_ACTIVE_SET_HPP

#include "utopia_Traits.hpp"

#include "utopia_BoxConstraints.hpp"
#include "utopia_DeviceView.hpp"

namespace utopia {
    template <class Vector>
    class ActiveSet {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using Layout = typename Traits<Vector>::Layout;

        virtual ~ActiveSet() = default;

        virtual bool determine(const BoxConstraints<Vector> &c, const Vector &x) {
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

    private:
        Vector indicator_;
        Scalar tol_{0};
        bool verbose_{false};
    };
}  // namespace utopia

#endif  // UTOPIA_ACTIVE_SET_HPP
