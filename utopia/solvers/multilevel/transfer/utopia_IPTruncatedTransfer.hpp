#ifndef UTOPIA_IP_MATRIX_TRUNCATED_TRANSFER_HPP
#define UTOPIA_IP_MATRIX_TRUNCATED_TRANSFER_HPP

#include "utopia_Transfer.hpp"
#include "utopia_TransferCommons.hpp"

#include <cassert>
#include <cmath>
#include <memory>

namespace utopia {
    /**
     * @brief      The class for transfer operators.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <class Matrix, class Vector>
    class IPTruncatedTransfer final : public MatrixTransfer<Matrix, Vector> {
    public:
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        using Scalar = typename utopia::Traits<Vector>::Scalar;

        IPTruncatedTransfer(const std::shared_ptr<Matrix> &I) {
            assert(I);
            I_ = I;
        }

        IPTruncatedTransfer(const std::shared_ptr<Matrix> &I, const std::shared_ptr<Matrix> &P) : I_(I), Pr_(P) {
            assert(I);
            assert(P);
        }

        ~IPTruncatedTransfer() override = default;

        /*=====================================================
                            initialization
        =====================================================*/

        /**
         * @brief      Initialization of interpolation operator.
         *
         * @param[in]  I_in  The interpolation.
         *
         */
        bool I_init(const std::shared_ptr<Matrix> &I_in) {
            assert(I_in);
            I_ = I_in;
            return true;
        }

        bool scale_transfer(const Matrix &scaling_mat) {
            assert(I_);
            *I_ = scaling_mat * (*I_);
            return true;
        }

        /*=====================================================
                                actions
        =====================================================*/
        /**
         * @brief      Interpolation of vector.
         *             \f$  x_{new} = I * x \f$
         *
         * @param[in]  x      The vector x.
         * @param      x_new  The interpoalted vector.
         *
         */
        bool interpolate(const Vector &x, Vector &x_new) const override {
            assert(I_);

            if (I_truncated_) {
                x_new = *I_truncated_ * x;
            } else {
                x_new = *I_ * x;
            }

            return true;
        }

        /**
         * @brief      Restriction of vector.
         *            \f$  x_{new} = R * x  \f$
         * @param[in]  x
         * @param      x_new
         *
         */
        bool restrict(const Vector &x, Vector &x_new) const override {
            assert(I_);

            if (I_truncated_) {
                x_new = transpose(*I_truncated_) * x;
            } else {
                x_new = transpose(*I_) * x;
            }
            return true;
        }

        /**
         * @brief      Restriction of vector based on booleans.
         * values have to be exact 1. and 0.
         *            \f$  x_{new} = R * x  \f$
         * @param[in]  x
         * @param      x_new
         *
         */
        bool boolean_restrict_or(const Vector &x, Vector &x_new) override {
            assert(I_);

            if (!I_truncated_) {
                Matrix R = transpose(*I_truncated_);
                utopia::boolean_restrict_or(R, x, x_new);
            } else {
                Matrix R = transpose(*I_);
                utopia::boolean_restrict_or(R, x, x_new);
            }

            return true;
        }

        /**
         * @brief      Restriction of matrix.
         *
         *             \f$  M_{new} = I^{T} * M  * I  \f$
         * @param[in]  M
         * @param      M_new
         *
         */
        bool restrict(const Matrix &M, Matrix &M_new) const override {
            assert(I_);

            if (I_truncated_) {
                M_new = utopia::ptap(M, *I_truncated_);
            } else {
                M_new = utopia::ptap(M, *I_);
            }

            return true;
        }

        /**
         * @brief      Initialization of projection down operator.
         *
         * @param[in]  P_in  The projection operator.
         *
         */
        bool P_init(const std::shared_ptr<Matrix> &P_in) {
            assert(P_in);
            Pr_ = P_in;
            return true;
        }

        /**
         * @brief      Projection of vector
         *            \f$  x_{new} = P * x  \f$
         * @param[in]  x
         * @param      x_new
         *
         */
        bool project_down(const Vector &x, Vector &x_new) const override {
            assert(Pr_);
            x_new = *Pr_ * x;
            return true;
        }

        bool project_down_positive_negative(const Vector &x_pos, const Vector &x_neg, Vector &x_new) override {
            if (empty(P_pos_)) {
                P_pos_ = *Pr_;
                chop_smaller_than(P_pos_, 1e-13);
            }

            if (empty(P_neg_)) {
                P_neg_ = (*Pr_);
                chop_greater_than(P_neg_, -1e-13);
            }

            x_new = (P_pos_ * x_pos) + (P_neg_ * x_neg);
            return true;
        }

        const Matrix &I() override {
            assert(I_);
            return *I_;
        }

        const Matrix &P() override {
            assert(Pr_);
            return *Pr_;
        }

        const Matrix &R() override {
            assert(false);
            Utopia::Abort("FIXME storing the restriction should not be a requirement!");
            return I();
        }

        Scalar interpolation_inf_norm() const override { return norm_infty(*I_); }

        Scalar projection_inf_norm() const override { return norm_infty(transpose(*I_)); }

        Scalar restriction_inf_norm() const override { return norm_infty(*Pr_); }

        void truncate_interpolation(const Vector &eq_active_flg) {
            // to speed up, we should check if constraint was changed between iterations (outside the method)
            if (!I_truncated_) {
                I_truncated_ = std::make_shared<Matrix>();
            }

            std::vector<SizeType> indices_eq_constraints_;
            {
                Read<Vector> r(eq_active_flg);

                Range range_w = range(eq_active_flg);
                for (SizeType i = range_w.begin(); i != range_w.end(); i++) {
                    if (eq_active_flg.get(i) == 1.0) {
                        indices_eq_constraints_.push_back(i);
                    }
                }
            }

            *I_truncated_ = *I_;
            set_zero_rows(*I_truncated_, indices_eq_constraints_, 0.0);
        }

        std::shared_ptr<Matrix> I_ptr() {
            assert(I_);
            return I_;
        }

        std::shared_ptr<Matrix> I_truncated_ptr() {
            assert(I_truncated_);
            return I_truncated_;
        }

    private:
        std::shared_ptr<Matrix> I_;
        std::shared_ptr<Matrix> Pr_;  // used only for nonlinear mutlilevel solvers
        Matrix P_pos_;                // used only for nonlinear mutlilevel solvers
        Matrix P_neg_;                // used only for nonlinear mutlilevel solvers

        std::shared_ptr<Matrix> I_truncated_;
    };

}  // namespace utopia

#endif  // UTOPIA_MATRIX_TRUNCATED_TRANSFER_HPP
