#ifndef UTOPIA_IPR_TRANSFER_HPP
#define UTOPIA_IPR_TRANSFER_HPP

#include "utopia_Readable.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Transfer.hpp"
#include "utopia_Writable.hpp"

#include "utopia_MatrixPtAPProduct.hpp"

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
    class IPRTransfer final : public MatrixTransfer<Matrix, Vector> {
    public:
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        using Scalar = typename utopia::Traits<Vector>::Scalar;

        IPRTransfer(const std::shared_ptr<Matrix> &I)
            :  // _I(I),
               // _R(transpose(I))
              I_norm_(0.0),
              R_norm_(0.0),
              P_norm_(0.0) {
            assert(I);

            _I = I;
            _R = std::make_shared<Matrix>(transpose(*I));
            _Pr = _R;
        }

        IPRTransfer(const std::shared_ptr<Matrix> &I, const std::shared_ptr<Matrix> &P)
            : _I(I), _R(std::make_shared<Matrix>(transpose(*I))), _Pr(P), I_norm_(0.0), R_norm_(0.0), P_norm_(0.0) {
            assert(I);
            assert(P);

            // std::cout<<"proper transfer down ... \n";
        }

        IPRTransfer(const std::shared_ptr<Matrix> &I,
                    const std::shared_ptr<Matrix> &R,
                    const std::shared_ptr<Matrix> &P)
            : _I(I), _R(R), _Pr(P), I_norm_(0.0), R_norm_(0.0), P_norm_(0.0) {
            assert(I);
            assert(R);
            assert(P);
        }

        ~IPRTransfer() override = default;

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
            assert(_R);
            x_new = *_R * x;
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

            Matrix R_boolean = *_R;

            R_boolean.transform_values(UTOPIA_LAMBDA(const Scalar value)->Scalar {
                return static_cast<Scalar>(device::abs(value) > off_diag_tol);
            });

            x_new = R_boolean * x;

            ReadAndWrite<Vector> rw_(x_new);

            auto r = range(x_new);
            for (auto i = r.begin(); i < r.end(); ++i) {
                if (x_new.get(i) > 1.) {
                    x_new.set(i, 1.);
                }
            }

            // THIS works for serial (use it once we have a is_parallel and is_serial
            // query) each_read(*_R, [&x, &x_new](const SizeType i, const SizeType j,
            // const Scalar value) {
            //     if(x.get(j) != 0. && std::abs(value) > off_diag_tol) {
            //         x_new.set(i, 1.);
            //     }
            // });

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

        const Matrix &I() override {
            assert(_I);
            return *_I;
        }

        const Matrix &R() override {
            assert(_R);
            return *_R;
        }

        const Matrix &P() override {
            assert(_Pr);
            return *_Pr;
        }

        Scalar interpolation_inf_norm() const override { return I_norm_ > 0 ? I_norm_ : norm_infty(*_I); }

        Scalar projection_inf_norm() const override {
            // return norm_infty(*_R);
            return P_norm_ > 0 ? P_norm_ : norm_infty(*_Pr);
        }

        Scalar restriction_inf_norm() const override {
            // return norm_infty(*_Pr);
            return R_norm_ > 0 ? R_norm_ : norm_infty(*_R);
        }

        void init_memory() override {
            I_norm_ = norm_infty(*_I);
            R_norm_ = norm_infty(*_R);
            P_norm_ = norm_infty(*_Pr);

            P_pos_ = *_Pr;
            chop_smaller_than(P_pos_, 1e-13);

            P_neg_ = (*_Pr);
            chop_greater_than(P_neg_, -1e-13);
        };

    private:
        std::shared_ptr<Matrix> _I, _R;  // _P;
        std::shared_ptr<Matrix> _Pr;
        Matrix P_pos_;
        Matrix P_neg_;
        Scalar I_norm_, R_norm_, P_norm_;
    };

}  // namespace utopia

#endif  // UTOPIA_MATRIX_TRANSFER_HPP
