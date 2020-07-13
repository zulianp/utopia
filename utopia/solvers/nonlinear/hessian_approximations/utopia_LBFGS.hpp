#ifndef UTOPIA_LBFGS_HPP
#define UTOPIA_LBFGS_HPP

#include "utopia_Core.hpp"
#include "utopia_HessianApproximation.hpp"

namespace utopia {
    enum LBFGSDampingTechnique { NOCEDAL = 0, POWEL = 1 };

    enum LBFGSScalingTechnique { NONE = 0, INITIAL = 1, ADAPTIVE = 2, FORBENIUS = 3 };

    // TODO:: speedup by preallocating just for H_inv, or H products...
    template <class Vector>
    class LBFGS final : public HessianApproximation<Vector> {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

    public:
        LBFGS(const SizeType &m)
            : m_(m),
              current_m_(0),
              theta_(1.0),
              gamma_(1.0),
              theta_min_(1.0),
              skip_update_treshold_(1e-12),
              damping_tech_(NOCEDAL),
              scaling_tech_(ADAPTIVE) {}

        ~LBFGS() override {
            Y_.clear();
            S_.clear();
            P_.clear();

            rho_.clear();
            yts_.clear();
            stp_.clear();

            SY_point_wise_.clear();
            YY_point_wise_.clear();
        }

        void initialize(const Vector &x_k, const Vector &g) override {
            HessianApproximation<Vector>::initialize(x_k, g);

            theta_ = 1.0;
            gamma_ = 1.0;

            current_m_ = 0;

            Y_.resize(m_);
            S_.resize(m_);
            P_.resize(m_);

            auto x_layout = layout(x_k);

            for (auto i = 0; i < m_; i++) {
                Y_[i].zeros(x_layout);
                S_[i].zeros(x_layout);
                P_[i].zeros(x_layout);
            }

            if (scaling_tech_ == LBFGSScalingTechnique::FORBENIUS) {
                SY_point_wise_.resize(m_);
                YY_point_wise_.resize(m_);

                for (auto i = 0; i < m_; i++) {
                    SY_point_wise_[i].zeros(x_layout);
                    YY_point_wise_[i].zeros(x_layout);
                }
            }

            D_.zeros(x_layout);
            D_inv_.zeros(x_layout);
            y_hat_.zeros(x_layout);
            help1_.zeros(x_layout);

            if (scaling_tech_ == LBFGSScalingTechnique::FORBENIUS) {
                SY_sum_.zeros(x_layout);
                YY_sum_.zeros(x_layout);
            }

            rho_.resize(m_);
            yts_.resize(m_);
            stp_.resize(m_);

            this->initialized(true);
        }

        void reset() override {
            theta_ = 1.0;
            gamma_ = 1.0;

            current_m_ = 0;

            // Y_.resize(m_);
            // S_.resize(m_);
            // P_.resize(m_);

            for (auto i = 0; i < m_; i++) {
                Y_[i].set(0.0);
                S_[i].set(0.0);
                P_[i].set(0.0);
            }

            if (scaling_tech_ == LBFGSScalingTechnique::FORBENIUS) {
                SY_point_wise_.resize(m_);
                YY_point_wise_.resize(m_);

                for (auto i = 0; i < m_; i++) {
                    SY_point_wise_[i].set(0.0);
                    YY_point_wise_[i].set(0.0);
                }
            }

            D_.set(0.0);
            D_inv_.set(0.0);
            y_hat_.set(0.0);
            help1_.set(0.0);

            if (scaling_tech_ == LBFGSScalingTechnique::FORBENIUS) {
                SY_sum_.set(0.0);
                YY_sum_.set(0.0);
            }

            // rho_.resize(m_);
            // yts_.resize(m_);
            // stp_.resize(m_);

            this->initialized(true);
        }

        inline LBFGS<Vector> *clone() const override { return new LBFGS<Vector>(*this); }

        bool update(const Vector &s, const Vector &y, const Vector & /*x*/, const Vector & /* g */) override {
            if (!this->initialized()) {
                utopia_error("BFGS::update: Initialization needs to be done before updating. \n");
                return false;
            }

            if (m_ == 0) {
                return true;
            }

            y_hat_ = y;

            // UTOPIA_NO_ALLOC_BEGIN("LBFGS1");
            bool skip_update = init_damping(y, s, y_hat_);
            // UTOPIA_NO_ALLOC_END();

            if (skip_update) {
                // UTOPIA_NO_ALLOC_BEGIN("LBFGS2");
                this->init_scaling_factors(y_hat_, s);
                // std::cout<<"updating..... \n";
                // UTOPIA_NO_ALLOC_END();
                return true;
            }

            Scalar denom = dot(y_hat_, s);

            if (current_m_ < m_) {
                // UTOPIA_NO_ALLOC_BEGIN("LBFGS3");
                Y_[current_m_] = y_hat_;
                S_[current_m_] = s;
                rho_[current_m_] = 1. / denom;

                yts_[current_m_] = denom;

                if (scaling_tech_ == LBFGSScalingTechnique::FORBENIUS) {
                    SY_point_wise_[current_m_] = e_mul(s, y_hat_);
                    YY_point_wise_[current_m_] = e_mul(y_hat_, y_hat_);
                }
                // UTOPIA_NO_ALLOC_END();
            } else {
                // UTOPIA_NO_ALLOC_BEGIN("LBFGS41");
                Y_[0] = y_hat_;
                S_[0] = s;
                rho_[0] = 1. / denom;
                yts_[0] = denom;
                // UTOPIA_NO_ALLOC_END();

                // UTOPIA_NO_ALLOC_BEGIN("LBFGS421");
                std::rotate(Y_.begin(), Y_.begin() + 1, Y_.end());
                std::rotate(S_.begin(), S_.begin() + 1, S_.end());
                std::rotate(rho_.begin(), rho_.begin() + 1, rho_.end());
                std::rotate(yts_.begin(), yts_.begin() + 1, yts_.end());
                // UTOPIA_NO_ALLOC_END();

                // UTOPIA_NO_ALLOC_BEGIN("LBFGS43");
                if (scaling_tech_ == LBFGSScalingTechnique::FORBENIUS) {
                    SY_point_wise_[0] = e_mul(s, y_hat_);
                    YY_point_wise_[0] = e_mul(y_hat_, y_hat_);

                    std::rotate(SY_point_wise_.begin(), SY_point_wise_.begin() + 1, SY_point_wise_.end());
                    std::rotate(YY_point_wise_.begin(), YY_point_wise_.begin() + 1, YY_point_wise_.end());
                }
                // UTOPIA_NO_ALLOC_END();
            }

            current_m_++;

            // updating the factors...
            // UTOPIA_NO_ALLOC_BEGIN("LBFGS5");
            this->init_scaling_factors(y_hat_, s);
            // UTOPIA_NO_ALLOC_END();

            // preallocation for forward product ...
            SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_;
            for (auto i = 0; i < current_memory_size; i++) {
                // UTOPIA_NO_ALLOC_BEGIN("LBFGS6");
                this->apply_H0(S_[i], P_[i]);
                // UTOPIA_NO_ALLOC_END();

                // UTOPIA_NO_ALLOC_BEGIN("LBFGS62");
                for (auto j = 0; j < i; j++) {
                    // Scalar mult1 = dot(S_[j], P_[i])/stp_[j];
                    // Scalar mult2 = dot(S_[i], Y_[j])/yts_[j];

                    Scalar mult1;
                    Scalar mult2;

                    dots(S_[j], P_[i], mult1, S_[i], Y_[j], mult2);
                    mult1 = mult1 / stp_[j];
                    mult2 = mult2 / yts_[j];

                    P_[i] = P_[i] - (mult1 * P_[j]) + (mult2 * Y_[j]);
                }
                // UTOPIA_NO_ALLOC_END();

                stp_[i] = dot(S_[i], P_[i]);
            }

            return true;
        }

        bool apply_Hinv(const Vector &g, Vector &q) override {
            if (!this->initialized()) {
                utopia_error("utopia::LBFGS::apply_Hinv:: missing initialization... \n");
            }

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

        bool apply_H(const Vector &v, Vector &result) override {
            return this - apply_H_new(v, result);

            // debugging
            // return this->apply_H_slow(v, result);
        }

        // version used for debugging purposes
        bool apply_H_slow(const Vector &v, Vector &result) {
            if (!this->initialized()) {
                utopia_error("utopia::LBFGS::apply_H:: missing initialization... \n");
            }

            std::vector<Vector> a, b;
            a.resize(m_);
            b.resize(m_);

            SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_;

            for (auto k = 0; k < current_memory_size; k++) {
                Scalar mult_b = 1. / std::sqrt(dot(Y_[k], S_[k]));
                b[k] = mult_b * Y_[k];

                this->apply_H0(S_[k], a[k]);

                for (auto i = 0; i < k; i++) {
                    Scalar bS = dot(b[i], S_[k]);
                    Scalar aS = dot(a[i], S_[k]);

                    a[k] += bS * b[i];
                    a[k] -= aS * a[i];
                }

                Scalar mult_a = 1. / std::sqrt(dot(S_[k], a[k]));
                a[k] = mult_a * a[k];
            }

            this->apply_H0(v, result);

            for (auto i = 0; i < current_memory_size; i++) {
                Scalar bv = dot(b[i], v);
                Scalar av = dot(a[i], v);

                result += (bv * b[i]);
                result -= (av * a[i]);
            }

            return true;
        }

        bool apply_H_new(const Vector &v, Vector &result) {
            if (!this->initialized()) {
                utopia_error("utopia::LBFGS::apply_H:: missing initialization... \n");
            }

            SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_;

            this->apply_H0(v, result);

            for (auto i = 0; i < current_memory_size; i++) {
                Scalar stz, ytx;
                dots(S_[i], result, stz, Y_[i], v, ytx);

                result = result + ((ytx / yts_[i]) * Y_[i]) - ((stz / stp_[i]) * P_[i]);
            }

            return true;
        }

        void memory_size(const SizeType &m) { m_ = m; }

        SizeType memory_size() const { return m_; }

        Scalar theta_min() const { return theta_min_; }

        LBFGSDampingTechnique damping_tech() const { return damping_tech_; }

        void damping_tech(const LBFGSDampingTechnique &damping) { damping_tech_ = damping; }

        LBFGSScalingTechnique scaling_tech() const { return scaling_tech_; }

        void scaling_tech(const LBFGSScalingTechnique &scaling) { scaling_tech_ = scaling; }

        void theta_min(const Scalar &theta_min) { theta_min_ = theta_min; }

        void read(Input &in) override {
            HessianApproximation<Vector>::read(in);
            in.get("memory_size", m_);
        }

        void print_usage(std::ostream &os) const override {
            HessianApproximation<Vector>::print_usage(os);
            this->print_param_usage(os, "memory_size", "int", "History/Memory size.", "5");
        }

        void H0_action(const std::function<void(const Vector &, Vector &)> operator_action) {
            H0_action_ = operator_action;
        }

        void H0_inv_action(const std::function<void(const Vector &, Vector &)> operator_action) {
            H0_inv_action_ = operator_action;
        }

    private:
        bool init_damping(const Vector &y, const Vector &s, Vector &y_hat) {
            if (damping_tech_ == LBFGSDampingTechnique::POWEL) {
                return init_damping_powel(y, s, y_hat);
            } else {
                return init_damping_nocedal(y, s, y_hat);
            }
        }

    public:  // nvcc requires any lambda to be public
        void init_scaling_factors(const Vector &y, const Vector &s) {
            if (scaling_tech_ == LBFGSScalingTechnique::NONE) {
                theta_ = 1.0;
                gamma_ = 1.0;

                D_.set(1.0);
                D_inv_.set(1.0);

            } else if (scaling_tech_ == LBFGSScalingTechnique::INITIAL) {
                if (current_m_ == 1) {
                    theta_ = dot(y, y) / dot(y, s);
                    theta_ = std::max(theta_min_, theta_);
                    gamma_ = 1. / theta_;

                    if (!std::isfinite(theta_) || !std::isfinite(gamma_)) {
                        theta_ = 1.0;
                        gamma_ = 1.0;
                    }

                    D_.set(theta_);
                    D_inv_.set(gamma_);
                }
            } else if (scaling_tech_ == LBFGSScalingTechnique::ADAPTIVE) {
                theta_ = dot(y, y) / dot(y, s);
                theta_ = std::max(theta_min_, theta_);
                gamma_ = 1. / theta_;

                if (!std::isfinite(theta_) || !std::isfinite(gamma_)) {
                    theta_ = 1.0;
                    gamma_ = 1.0;
                }

                D_.set(theta_);
                D_inv_.set(gamma_);
            }
            // FORBENIUS
            else {
                SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_;

                SY_sum_.set(0.0);
                YY_sum_.set(0.0);

                for (auto i = 0; i < current_memory_size; i++) {
                    SY_sum_ += SY_point_wise_[i];
                    YY_sum_ += YY_point_wise_[i];
                }

                Scalar theta_working = dot(y, y) / dot(y, s);
                theta_working = std::max(theta_min_, theta_working);

                if (!std::isfinite(theta_working)) {
                    theta_working = 1.0;
                }

                {
                    auto d_YY_sum = const_local_view_device(YY_sum_);
                    auto d_SY_sum = const_local_view_device(SY_sum_);

                    auto D_view = local_view_device(D_);

                    parallel_for(
                        local_range_device(D_), UTOPIA_LAMBDA(const SizeType i) {
                            Scalar yy = d_YY_sum.get(i);
                            Scalar sy = d_SY_sum.get(i);

                            if (!device::isfinite(yy) || !device::isfinite(sy)) {
                                D_view.set(i, theta_working);
                            }

                            Scalar dii = (yy / sy);

                            if (dii > 10e-10 && device::isfinite(dii)) {
                                D_view.set(i, dii);
                            } else {
                                D_view.set(i, theta_working);
                            }
                        });
                }

                {
                    auto D_view = const_local_view_device(D_);
                    auto D_inv_view = local_view_device(D_inv_);

                    parallel_for(
                        local_range_device(D_), UTOPIA_LAMBDA(const SizeType i) {
                            Scalar dii = D_view.get(i);

                            if (device::isfinite(dii) && dii != 0.0) {
                                D_inv_view.set(i, 1. / dii);
                            } else {
                                D_inv_view.set(i, theta_working);
                            }
                        });
                }

            }  // end FORM norm

        }  // init_scaling_factors()

    private:
        bool init_damping_nocedal(const Vector &y, const Vector &s, Vector &y_hat) {
            bool skip_update = false;
            if (dot(y, s) < skip_update_treshold_) {
                // if(mpi_world_rank()==0){
                //     utopia_warning("L-BFGS-B: Curvature condition not satified. Skipping update. \n");
                // }

                skip_update = true;
            }

            y_hat = y;

            return skip_update;
        }

        bool init_damping_powel(const Vector &y, const Vector &s, Vector &y_hat) {
            // FIXME remove temporaries
            // Vector Bs = 0.0*s;

            // help1_ = 0.0*s;
            help1_.set(0.0);
            this->apply_H(s, help1_);

            // Scalar sy   = dot(y,s);
            // Scalar sBs  = dot(s, Bs);

            Scalar sy, sBs;
            dots(y, s, sy, s, help1_, sBs);

            Scalar psi;

            if (sy >= 0.2 * sBs) {
                psi = 1.0;
            } else {
                psi = (0.8 * sBs) / (sBs - sy);
            }

            y_hat = (psi * y) + ((1.0 - psi) * help1_);

            return false;
        }

        void apply_H0(const Vector &v, Vector &result) {
            if (H0_action_ && current_m_ == 0) {
                // UTOPIA_NO_ALLOC_BEGIN("LBGFS:A1");
                H0_action_(v, result);
                // UTOPIA_NO_ALLOC_END();
            } else {
                // why is this different ???
                if (current_m_ == 0) {
                    // UTOPIA_NO_ALLOC_BEGIN("LBGFS:A2");
                    result = theta_ * v;
                    // UTOPIA_NO_ALLOC_END();
                } else {
                    // UTOPIA_NO_ALLOC_BEGIN("LBGFS:A22");
                    result = e_mul(D_, v);
                    // UTOPIA_NO_ALLOC_END();
                }
            }
        }

        void apply_H0_inv(const Vector &v, Vector &result) {
            if (H0_inv_action_ && current_m_ == 0) {
                H0_inv_action_(v, result);
            } else {
                if (current_m_ == 0) {
                    result = gamma_ * v;
                } else {
                    result = e_mul(D_inv_, v);
                }
            }
        }

        Scalar compute_uHu_dot(const Vector &u) override {
            this->apply_H(u, help1_);
            return dot(u, help1_);
        }

        Scalar compute_uHinvv_dot(const Vector &u, const Vector &v) override {
            this->apply_Hinv(v, help1_);
            return dot(u, help1_);
        }

        Scalar compute_uHv_dot(const Vector &u, const Vector &v) override {
            this->apply_H(v, help1_);
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
        std::vector<Vector> P_;

        std::vector<Scalar> rho_;
        std::vector<Scalar> yts_;
        std::vector<Scalar> stp_;

        LBFGSDampingTechnique damping_tech_;
        LBFGSScalingTechnique scaling_tech_;

        Vector D_;
        Vector D_inv_;
        Vector y_hat_;
        Vector help1_;

        Vector SY_sum_;
        Vector YY_sum_;

        std::vector<Vector> SY_point_wise_;
        std::vector<Vector> YY_point_wise_;

        std::function<void(const Vector &, Vector &)> H0_action_;
        std::function<void(const Vector &, Vector &)> H0_inv_action_;
    };

}  // namespace utopia

#endif  // UTOPIA_LBFGS_HPP
