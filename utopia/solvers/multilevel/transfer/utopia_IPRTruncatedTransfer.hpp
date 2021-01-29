#ifndef UTOPIA_IPR_MATRIX_TRUNCATED_TRANSFER_HPP
#define UTOPIA_IPR_MATRIX_TRUNCATED_TRANSFER_HPP

#include "utopia_Transfer.hpp"

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
    class IPRTruncatedTransfer final : public MatrixTransfer<Matrix, Vector> {
    public:
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        using Scalar = typename utopia::Traits<Vector>::Scalar;

        IPRTruncatedTransfer(const std::shared_ptr<Matrix> &I) {
            assert(I);

            _I = I;
            _R = std::make_shared<Matrix>(transpose(*I));
            // Pr_ = _R;

            R_truncated_ = std::make_shared<Matrix>(*_R);
            I_truncated_ = std::make_shared<Matrix>(*_I);
        }

        IPRTruncatedTransfer(const std::shared_ptr<Matrix> &I, const std::shared_ptr<Matrix> &P)
            : _I(I), _R(std::make_shared<Matrix>(transpose(*I))), Pr_(P) {
            assert(I);
            assert(P);

            R_truncated_ = std::make_shared<Matrix>(*_R);
            I_truncated_ = std::make_shared<Matrix>(*_I);

            // std::cout<<"proper transfer down ... \n";
        }

        IPRTruncatedTransfer(const std::shared_ptr<Matrix> &I,
                             const std::shared_ptr<Matrix> &R,
                             const std::shared_ptr<Matrix> &P)
            : _I(I), _R(R), Pr_(P) {
            assert(I);
            assert(R);
            assert(P);
        }

        ~IPRTruncatedTransfer() override = default;

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
            _I = I_in;
            return true;
        }

        /**
         * @brief      Initialization of restriction operator.
         *
         * @param[in]  R_in  The restriction.
         *
         */
        bool R_init(const std::shared_ptr<Matrix> &R_in) {
            assert(R_in);
            _R = R_in;
            return true;
        }

        /**
         * @brief      Initialization of interpolation & restriction operators.
         *
         * @param[in]  I_in  The projection.
         * @param[in]  R_in  The restriction.
         *
         * @return
         */
        bool IR_init(const std::shared_ptr<Matrix> &I_in, const std::shared_ptr<Matrix> &R_in) {
            assert(I_in);
            assert(R_in);

            _I = I_in;
            _R = R_in;

            return true;
        }

        bool scale_transfer(const Matrix &scaling_mat) {
            assert(_I);
            assert(_R);

            *_I = scaling_mat * (*_I);
            *_R = transpose(*_I);

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
            assert(I_truncated_);
            x_new = *I_truncated_ * x;
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
            assert(R_truncated_);
            x_new = *R_truncated_ * x;
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
            static const Scalar off_diag_tol = std::numeric_limits<Scalar>::epsilon() * 1e6;
            if (!R_truncated_) {
                *R_truncated_ = *_R;
            }

            Matrix R_boolean = *R_truncated_;

            R_boolean.transform_values(UTOPIA_LAMBDA(const Scalar &value)->Scalar {
                if (device::abs(value) > off_diag_tol) {
                    return 1.;
                } else {
                    return 0.;
                }
            });

            x_new = R_boolean * x;

            x_new.transform_values(UTOPIA_LAMBDA(const Scalar &value)->Scalar {
                if (value > 1.) {
                    return 1.0;
                } else {
                    return 0.0;
                }
            });

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
            assert(I_truncated_);
            M_new = utopia::ptap(M, *I_truncated_);
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
            assert(_I);
            return *_I;
        }

        const Matrix &R() override {
            assert(_R);
            return *_R;
        }
        const Matrix &P() override {
            assert(Pr_);
            return *Pr_;
        }

        Scalar interpolation_inf_norm() const override { return norm_infty(*_I); }

        Scalar projection_inf_norm() const override { return norm_infty(*_R); }

        Scalar restriction_inf_norm() const override { return norm_infty(*Pr_); }

        void truncate_interpolation(const Vector &eq_active_flg) {
            // to speed up, we should check if constraint was changed between iterations
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

            *I_truncated_ = *_I;
            set_zero_rows(*I_truncated_, indices_eq_constraints_, 0.0);
            *R_truncated_ = transpose(*I_truncated_);
        }

        void detect_zero_rows_on_coarser_level(const Vector &eq_active_flg, Vector &zero_rows) const {
            // Propagate flags to coarser level
            zero_rows = *R_truncated_ * eq_active_flg;

            auto zero_rows_view = local_view_device(zero_rows);

            parallel_for(
                local_range_device(zero_rows), UTOPIA_LAMBDA(const SizeType i) {
                    auto val = zero_rows_view.get(i);
                    if (device::abs(val) < 1e-16) {
                        // numerical zero
                        val = 1.0;
                    } else {
                        val = 0.0;
                    }

                    zero_rows_view.set(val);
                });
        }

    private:
        std::shared_ptr<Matrix> _I, _R;

        std::shared_ptr<Matrix> Pr_;  // used only for nonlinear mutlilevel solvers
        Matrix P_pos_;                // used only for nonlinear mutlilevel solvers
        Matrix P_neg_;                // used only for nonlinear mutlilevel solvers

        std::shared_ptr<Matrix> I_truncated_;
        std::shared_ptr<Matrix> R_truncated_;
    };

}  // namespace utopia

#endif  // UTOPIA_MATRIX_TRUNCATED_TRANSFER_HPP
