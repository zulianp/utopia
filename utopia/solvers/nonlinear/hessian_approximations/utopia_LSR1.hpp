#ifndef UTOPIA_LSR1_HPP
#define UTOPIA_LSR1_HPP

#include "utopia_Core.hpp"
#include "utopia_HessianApproximation.hpp"

namespace utopia {
    template <class Vector>
    class LSR1 final : public HessianApproximation<Vector> {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

    public:
        LSR1(const SizeType &m) : m_(m), current_m_(0), theta_(1.0), gamma_(1.0) {}

        ~LSR1() override {
            Y_.clear();
            S_.clear();

            p_.clear();
            p_inv_.clear();
        }

        void initialize(const Vector &x_k, const Vector &g) override {
            HessianApproximation<Vector>::initialize(x_k, g);

            theta_ = 1.0;
            gamma_ = 1.0;
            current_m_ = 0;

            Y_.resize(m_);
            S_.resize(m_);

            p_.resize(m_);
            p_inv_.resize(m_);

            auto x_layout = layout(x_k);

            for (auto i = 0; i < m_; i++) {
                Y_[i].zeros(x_layout);
                S_[i].zeros(x_layout);
                p_[i].zeros(x_layout);
                p_inv_[i].zeros(x_layout);
            }

            this->initialized(true);
        }

        void reset() override {
            // Y_.clear();
            // S_.clear();

            // p_.clear();
            // p_inv_.clear();

            // Vector x, g;
            // this->initialize(x, g);

            theta_ = 1.0;
            gamma_ = 1.0;
            current_m_ = 0;

            Y_.resize(m_);
            S_.resize(m_);

            p_.resize(m_);
            p_inv_.resize(m_);

            for (auto i = 0; i < m_; i++) {
                Y_[i].set(0.0);
                S_[i].set(0.0);
                p_[i].set(0.0);
                p_inv_[i].set(0.0);
            }

            this->initialized(true);
        }

        inline LSR1<Vector> *clone() const override { return new LSR1<Vector>(*this); }

        bool update(const Vector &s, const Vector &y, const Vector & /*x*/, const Vector & /* g */) override {
            if (!this->initialized()) {
                utopia_error("BFGS::update: Initialization needs to be done before updating. \n");
                return false;
            } else if (m_ == 0) {
                return true;
            }

            diff_ = y - s;
            // Scalar nom = std::abs(dot(s, diff_));
            // Scalar denom = norm2(s) * norm2(diff);

            Scalar nom, denom, denom2;
            dots(s, diff_, nom, s, s, denom, diff_, diff_, denom2);

            nom = std::abs(nom);
            denom = denom * denom2;

            if (nom / denom < this->num_tol() || !std::isfinite(denom) || !std::isfinite(nom)) {
                // if(mpi_world_rank()==0 ){
                //     utopia_warning("L-SR1: Curvature condition not satified. Skipping
                //     update. \n");
                // }

                return false;
            }

            if (current_m_ < m_) {
                Y_[current_m_] = y;
                S_[current_m_] = s;
            } else {
                Y_[0] = y;
                S_[0] = s;

                std::rotate(Y_.begin(), Y_.begin() + 1, Y_.end());
                std::rotate(S_.begin(), S_.begin() + 1, S_.end());
            }

            current_m_++;

            precompute_p();
            precompute_p_inv();

            return true;
        }

        bool apply_Hinv(const Vector &v, Vector &result) override {
            if (!this->initialized()) {
                utopia_error("utopia::LSR1::apply_Hinv:: missing initialization... \n");
            }

            SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_;
            result = gamma_ * v;

            for (auto i = 0; i < current_memory_size; i++) {
                Scalar scaling_factor = dot(p_inv_[i], v) / dot(p_inv_[i], Y_[i]);
                result += scaling_factor * p_inv_[i];
            }

            if (has_nan_or_inf(result)) result = gamma_ * v;

            return true;
        }

        bool apply_H(const Vector &v, Vector &result) override {
            if (!this->initialized()) {
                utopia_error("utopia::LSR1::apply_Hinv:: missing initialization... \n");
            }

            SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_;
            result = theta_ * v;

            for (auto i = 0; i < current_memory_size; i++) {
                Scalar scaling_factor = dot(p_[i], v) / dot(p_[i], S_[i]);
                result += scaling_factor * p_[i];
            }

            if (has_nan_or_inf(result)) result = theta_ * v;

            return true;
        }

        void memory_size(const SizeType &m) { m_ = m; }

        SizeType memory_size() const override { return m_; }

        void read(Input &in) override {
            HessianApproximation<Vector>::read(in);
            in.get("memory_size", m_);
        }

        void print_usage(std::ostream &os) const override {
            HessianApproximation<Vector>::print_usage(os);
            this->print_param_usage(os, "memory_size", "int", "History/Memory size.", "5");
        }

    private:
        void precompute_p() {
            if (!this->initialized()) {
                utopia_error("utopia::LSR1::apply_Hinv:: missing initialization... \n");
            }

            SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_;

            for (auto i = 0; i < current_memory_size; i++) {
                p_[i] = Y_[i] - theta_ * S_[i];
                for (auto j = 0; j < i; j++) {
                    Scalar scaling_factor = dot(p_[j], S_[i]) / dot(p_[j], S_[j]);
                    p_[i] -= scaling_factor * p_[j];
                }
            }
        }

        void precompute_p_inv() {
            if (!this->initialized()) {
                utopia_error("utopia::LSR1::apply_Hinv:: missing initialization... \n");
            }

            SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_;

            for (auto i = 0; i < current_memory_size; i++) {
                p_inv_[i] = S_[i] - gamma_ * Y_[i];
                for (auto j = 0; j < i; j++) {
                    Scalar scaling_factor = dot(p_inv_[j], Y_[i]) / dot(p_inv_[j], Y_[j]);
                    p_inv_[i] -= scaling_factor * p_inv_[j];
                }
            }
        }

        Scalar compute_uHu_dot(const Vector &u) override {
            this->apply_H(u, diff_);
            return dot(u, diff_);
        }

        Scalar compute_uHinvv_dot(const Vector &u, const Vector &v) override {
            this->apply_Hinv(v, diff_);
            return dot(u, diff_);
        }

        Scalar compute_uHv_dot(const Vector &u, const Vector &v) override {
            this->apply_H(v, diff_);
            return dot(u, diff_);
        }

    private:
        SizeType m_;          // memory size
        SizeType current_m_;  // current amount of vectors in the memory

        Scalar theta_;
        Scalar gamma_;

        std::vector<Vector> Y_;
        std::vector<Vector> S_;

        std::vector<Vector> p_;
        std::vector<Vector> p_inv_;

        Vector diff_;
    };

}  // namespace utopia

#endif  // UTOPIA_LBFGS_HPP
