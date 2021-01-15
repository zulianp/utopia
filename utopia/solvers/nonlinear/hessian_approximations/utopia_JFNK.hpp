#ifndef UTOPIA_JFNK_HPP
#define UTOPIA_JFNK_HPP

#include "utopia_Core.hpp"
#include "utopia_HessianApproximation.hpp"

namespace utopia {

    template <class Vector>
    class JFNK final : public HessianApproximation<Vector> {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

    public:
        JFNK(const FunctionBase<Vector> &fun) : eps_(1e-14), fun_(fun) {}

        void initialize(const Vector &x_k, const Vector &g) override {
            HessianApproximation<Vector>::initialize(x_k, g);

            x_k_ = x_k;
            g_ = g;

            auto vec_lo = layout(g);

            grad_pertubed_.values(vec_lo, 0.0);
            ones_.values(vec_lo, 1.0);
            x_p_.values(vec_lo, 0.0);

            this->initialized(true);
        }

        void reset() override {
            Vector x, g;
            this->initialize(x, g);
        }

        inline JFNK<Vector> *clone() const override { return new JFNK<Vector>(*this); }

        bool update(const Vector & /*s*/, const Vector & /*y*/, const Vector &x, const Vector &g) override {
            x_k_ = x;
            g_ = g;

            return true;
        }

        bool apply_Hinv(const Vector & /*g*/, Vector & /*q*/) override {
            utopia_error("utopia::JFNK::apply_Hinv:: not supported... \n");
            return false;
        }

        // TODO:: make more memory efficient
        bool apply_H(const Vector &v, Vector &result) override {
            x_p_ = std::sqrt(eps_) * (ones_ + x_k_);
            Scalar sum_a = sum(x_p_);

            SizeType n_glob = size(result).get(0);

            Scalar per_ = 0.0;

            if (norm2(v) > eps_) {
                per_ = 1. / (n_glob * Scalar(norm2(v)));
                per_ *= sum_a;
            } else {
                per_ = sum_a / n_glob;
            }

            x_p_ = x_k_ + (per_ * v);
            fun_.gradient(x_p_, grad_pertubed_);

            result = (1. / per_) * (grad_pertubed_ - g_);

            return true;
        }

        void read(Input &in) override { HessianApproximation<Vector>::read(in); }

        void print_usage(std::ostream &os) const override { HessianApproximation<Vector>::print_usage(os); }

    private:
        Scalar eps_;
        const FunctionBase<Vector> &fun_;
        Vector x_k_, g_, ones_, grad_pertubed_, x_p_;
    };

}  // namespace utopia

#endif  // UTOPIA_JFNK_HPP
