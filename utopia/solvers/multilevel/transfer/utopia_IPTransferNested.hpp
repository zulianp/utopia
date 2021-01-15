#ifndef UTOPIA_MATRIX_TRANSFER_NESTED_HPP
#define UTOPIA_MATRIX_TRANSFER_NESTED_HPP

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
    class IPTransferNested final : public MatrixTransfer<Matrix, Vector> {
    public:
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        using Scalar = typename utopia::Traits<Vector>::Scalar;

        IPTransferNested(const std::shared_ptr<Matrix> &I, const std::shared_ptr<Matrix> &P)
            : _I(I),
              // _R(std::make_shared<Matrix>(transpose(*I))),
              _Pr(P),
              I_norm_(1.0),
              R_norm_(1.0),
              P_norm_(1.0) {
            assert(I);
            assert(P);
        }

        ~IPTransferNested() override = default;

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
            return true;
        }

        // const Matrix &R() const
        // {
        //     return *_R;
        // }

        /**
         * @brief      Restriction of vector based on booleans.
         * values have to be exact 1. and 0.
         *            \f$  x_{new} = R * x  \f$
         * @param[in]  x
         * @param      x_new
         *
         */
        bool boolean_restrict_or(const Vector & /*x*/, Vector & /*x_new*/) override {
            assert(false && "implement me");
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
            assert(_I);
            M_new = utopia::ptap(M, *_I);
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
            _Pr = P_in;
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
            assert(_Pr);
            x_new = *_Pr * x;
            return true;
        }

        bool project_down_positive_negative(const Vector &x_pos, const Vector & /*x_neg*/, Vector &x_new) override {
            x_new = *_Pr * x_pos;
            return true;
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

        void handle_equality_constraints(const Vector &is_constrained) override {
            assert(_I);
            assert(!empty(*_I));
            assert(!empty(is_constrained));
            assert(size(is_constrained).get(0) == size(*_I).get(0));

            if (empty(is_constrained) || empty(*_I)) return;

            set_zero_rows(*_I, is_constrained, 0.0);
        }

        Scalar interpolation_inf_norm() const override { return I_norm_ > 0 ? I_norm_ : norm_infty(*_I); }

        Scalar projection_inf_norm() const override { return P_norm_ > 0 ? P_norm_ : norm_infty(*_Pr); }

        Scalar restriction_inf_norm() const override { return R_norm_ > 0 ? R_norm_ : norm_infty(transpose(*_I)); }

        void init_memory() override {
            I_norm_ = norm_infty(*_I);
            R_norm_ = norm_infty(transpose(*_I));
            P_norm_ = norm_infty(*_Pr);
        };

    private:
        std::shared_ptr<Matrix> _I;  // _R;
        std::shared_ptr<Matrix> _Pr;
        std::shared_ptr<Matrix> _R;
        Scalar I_norm_, R_norm_, P_norm_;
    };

}  // namespace utopia

#endif  // UTOPIA_MATRIX_TRANSFER_NESTED_HPP
