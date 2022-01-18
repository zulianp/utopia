#ifndef UTOPIA_LBFGS_INV_HPP
#define UTOPIA_LBFGS_INV_HPP

#include "utopia_Core.hpp"
#include "utopia_HessianApproximation.hpp"

namespace utopia {

    // Class considers computations only for inverse approximation
    // If you need also actuall Hessian, take LBFGS class
    template <class Vector>
    class LBFGSInv final : public HessianApproximation<Vector> {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

    public:
        LBFGSInv(const SizeType &m)
            : m_(m), current_m_(0), theta_(1.0), gamma_(1.0), theta_min_(1.0), skip_update_treshold_(1e-15) {}

        ~LBFGSInv() override {
            Y_.clear();
            S_.clear();
            rho_.clear();
        }

        SizeType current_m() { return current_m_; }

        bool is_approx_fully_built() {
            if (current_m_ < m_)
                return false;
            else
                return true;
        }

        void initialize(const Vector &x_k, const Vector &g) override {
            HessianApproximation<Vector>::initialize(x_k, g);

            theta_ = 1.0;
            gamma_ = 1.0;

            current_m_ = 0;

            Y_.resize(m_);
            S_.resize(m_);

            auto x_layout = layout(x_k);

            for (auto i = 0; i < m_; i++) {
                Y_[i].zeros(x_layout);
                S_[i].zeros(x_layout);
            }

            help1_.zeros(x_layout);

            rho_.resize(m_);
            this->initialized(true);
        }

        void reset() override {
            theta_ = 1.0;
            gamma_ = 1.0;

            current_m_ = 0;

            for (auto i = 0; i < m_; i++) {
                Y_[i].set(0.0);
                S_[i].set(0.0);
            }

            help1_.set(0.0);

            this->initialized(true);
        }

        inline LBFGSInv<Vector> *clone() const override { return new LBFGSInv<Vector>(*this); }

        bool update(const Vector &s, const Vector &y, const Vector & /*x*/, const Vector & /* g */) override {
            if (!this->initialized()) {
                utopia_error("BFGS::update: Initialization needs to be done before updating. \n");
                return false;
            }

            if (m_ == 0) {
                return true;
            }

            bool skip_update = init_damping(y, s);

            if (skip_update) {
                // utopia::out() << "update skipped! \n";
                this->init_scaling_factors(y, s);
                return true;
            }

            Scalar denom = dot(y, s);

            if (current_m_ < m_) {
                Y_[current_m_] = y;
                S_[current_m_] = s;
                rho_[current_m_] = 1. / denom;

            } else {
                Y_[0] = y;
                S_[0] = s;
                rho_[0] = 1. / denom;

                std::rotate(Y_.begin(), Y_.begin() + 1, Y_.end());
                std::rotate(S_.begin(), S_.begin() + 1, S_.end());
                std::rotate(rho_.begin(), rho_.begin() + 1, rho_.end());
            }

            current_m_++;

            // updating the factors...
            this->init_scaling_factors(y, s);

            return true;
        }

        bool replace_at_update_inv(const SizeType &index, const Vector &s, const Vector &y) override {
            if (!this->initialized()) {
                utopia_error("L-BFGS::update: Initialization needs to be done before updating. \n");
                return false;
            }

            if (index > m_ or index > current_m_) {
                utopia_error("L-BFGS::update: inndex does not exist  \n");
                return false;
            }

            bool skip_update = init_damping(y, s);

            if (skip_update) {
                this->init_scaling_factors(y, s);
                return true;
            }

            Scalar denom = dot(y, s);

            // replacing stuff
            Y_[index] = y;
            S_[index] = s;
            rho_[index] = 1. / denom;

            // updating the factors...
            this->init_scaling_factors(y, s);

            return true;
        }

        bool apply_Hinv(const Vector &g, Vector &q) override {
            if (!this->initialized()) {
                utopia_error("utopia::LBFGS::apply_Hinv:: missing initialization... \n");
            }

            // std::cout << "current_m_: " << current_m_ << "  \n";

            SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_;
            std::vector<Scalar> alpha_inv(current_memory_size);
            q = g;

            for (auto i = current_memory_size - 1; i >= 0; i--) {
                alpha_inv[i] = rho_[i] * dot(S_[i], q);
                q -= alpha_inv[i] * Y_[i];
            }

            help1_ = q;
            this->apply_H0_inv(help1_, q);

            for (auto i = 0; i < current_memory_size; i++) {
                Scalar betta_inv = rho_[i] * dot(Y_[i], q);
                q += (alpha_inv[i] - betta_inv) * S_[i];
            }

            if (has_nan_or_inf(q)) {
                this->apply_H0_inv(g, q);
            }

            return true;
        }

        bool apply_H(const Vector & /*v*/, Vector & /*result*/) override {
            utopia_error("L-BFGS_inv - apply_H not implemented... \n");
            return false;
        }

        void memory_size(const SizeType &m) { m_ = m; }

        SizeType memory_size() const override { return m_; }

        Scalar theta_min() const { return theta_min_; }

        void theta_min(const Scalar &theta_min) { theta_min_ = theta_min; }

        void read(Input &in) override {
            HessianApproximation<Vector>::read(in);
            in.get("memory_size", m_);
        }

        void print_usage(std::ostream &os) const override {
            HessianApproximation<Vector>::print_usage(os);
            this->print_param_usage(os, "memory_size", "int", "History/Memory size.", "5");
        }

        void H0_inv_action(const std::function<void(const Vector &, Vector &)> operator_action) {
            H0_inv_action_ = operator_action;
        }

    private:
        bool init_damping(const Vector &y, const Vector &s) { return init_damping_nocedal(y, s); }

    public:  // nvcc requires any lambda to be public
        void init_scaling_factors(const Vector &y, const Vector &s) {
            theta_ = dot(y, y) / dot(y, s);
            theta_ = std::max(theta_min_, theta_);
            gamma_ = 1. / theta_;

            if (!std::isfinite(theta_) || !std::isfinite(gamma_)) {
                gamma_ = 1.0;
            }

        }  // init_scaling_factors()

    private:
        bool init_damping_nocedal(const Vector &y, const Vector &s) {
            bool skip_update = false;

            // std::cout << "dot(y, s): " << dot(y, s) << "  \n";

            if (dot(y, s) < skip_update_treshold_) {
                // if(mpi_world_rank()==0){
                //     utopia_warning("L-BFGS-B: Curvature condition not satified.
                //     Skipping update. \n");
                // }

                skip_update = true;
            }
            return skip_update;
        }

        void apply_H0_inv(const Vector &v, Vector &result) {
            if (H0_inv_action_ && current_m_ == 0) {
                H0_inv_action_(v, result);
            } else {
                result = gamma_ * v;
            }
        }

        Scalar compute_uHinvv_dot(const Vector &u, const Vector &v) override {
            this->apply_Hinv(v, help1_);
            return dot(u, help1_);
        }

    private:
        SizeType m_;          // memory size
        SizeType current_m_;  // current amount of vectors in the memory

        Scalar theta_;
        Scalar gamma_;
        Scalar theta_min_;
        Scalar skip_update_treshold_;

        std::vector<Vector> Y_;
        std::vector<Vector> S_;

        std::vector<Scalar> rho_;
        Vector help1_;

        std::function<void(const Vector &, Vector &)> H0_inv_action_;
    };

}  // namespace utopia

#endif  // UTOPIA_LBFGS_INV_HPP
