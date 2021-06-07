#ifndef UTOPIA_IP_TRANSFER_HPP
#define UTOPIA_IP_TRANSFER_HPP

#include "utopia_Temp.hpp"
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
    class IPTransfer final : public MatrixTransfer<Matrix, Vector> {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

    public:
        IPTransfer(const std::shared_ptr<Matrix> &I, const Scalar restrict_factor = 1.)
            : restrict_factor_(restrict_factor) {
            assert(I);
            _I = I;
        }

        IPTransfer(const std::shared_ptr<Matrix> &I, const std::shared_ptr<Matrix> &P)
            : _I(I), _Pr(P), restrict_factor_(1.) {
            assert(I);
            assert(P);
        }

        ~IPTransfer() override = default;

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
            assert(_I);
            x_new = *_I * x;
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
            assert(_I);
            x_new = transpose(*_I) * x;

            if (restrict_factor_ != 1.0) {
                x_new *= restrict_factor_;
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
            Matrix R = transpose(*_I);
            return utopia::boolean_restrict_or(R, x, x_new);
        }

        void handle_equality_constraints(const Vector &is_constrained) override {
            assert(_I);
            assert(!empty(*_I));
            assert(!empty(is_constrained));
            assert(size(is_constrained).get(0) == size(*_I).get(0));

            if (empty(is_constrained) || empty(*_I)) return;

            set_zero_rows(*_I, is_constrained, 0.0);
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
            assert(_I);
            M_new = utopia::ptap(M, *_I);
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
            assert(_Pr || _I);

            if (!_Pr) {
                return restrict(x, x_new);
            } else {
                x_new = *_Pr * x;
            }
            return true;
        }

        bool project_down_positive_negative(const Vector &x_pos, const Vector &x_neg, Vector &x_new) override {
            if (empty(P_pos_)) {
                P_pos_ = *_Pr;
                chop_smaller_than(P_pos_, 1e-13);
            }

            if (empty(P_neg_)) {
                P_neg_ = (*_Pr);
                chop_greater_than(P_neg_, -1e-13);
            }

            x_new = (P_pos_ * x_pos) + (P_neg_ * x_neg);
            return true;
        }

        Scalar interpolation_inf_norm() const override { return norm_infty(*_I); }

        Scalar projection_inf_norm() const override {
            m_utopia_warning_once("projection_inf_norm not implemented properly");
            return norm_infty(transpose(*_I));
        }

        Scalar restriction_inf_norm() const override { return norm_infty(*_Pr); }

        std::shared_ptr<Matrix> I_ptr() {
            assert(_I);
            return _I;
        }

        const Matrix &I() override {
            assert(_I);
            return *_I;
        }

        const Matrix &P() override {
            assert(_Pr);
            return *_Pr;
        }

        const Matrix &R() override {
            assert(_I);
            _R = std::make_shared<Matrix>(transpose(*_I));
            return *_R;
        }

    protected:
        std::shared_ptr<Matrix> _I;
        std::shared_ptr<Matrix> _Pr;
        std::shared_ptr<Matrix> _R;

        Scalar restrict_factor_;

        Matrix P_pos_;
        Matrix P_neg_;
    };

}  // namespace utopia

#endif  // UTOPIA_IP_TRANSFER_HPP
