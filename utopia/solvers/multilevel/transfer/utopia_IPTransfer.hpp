#ifndef UTOPIA_IP_TRANSFER_HPP
#define UTOPIA_IP_TRANSFER_HPP

#include "utopia_Transfer.hpp"
#include "utopia_Temp.hpp"

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
        template<class Matrix, class Vector>
        class IPTransfer final : public Transfer<Matrix, Vector>
        {
            typedef UTOPIA_SCALAR(Vector)    Scalar;
            typedef UTOPIA_SIZE_TYPE(Vector) SizeType;


        public:

        IPTransfer(const std::shared_ptr<Matrix> &I, const Scalar restrict_factor = 1.)
        : restrict_factor_(restrict_factor)
        {
            assert(I);
            _I = I;
        }

        IPTransfer(const std::shared_ptr<Matrix> &I, const std::shared_ptr<Matrix> &P):
                _I(I),
                _Pr(P),
                restrict_factor_(1.)
        {
            assert(I);
            assert(P);
        }


         ~IPTransfer(){}


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
         bool interpolate(const Vector &x, Vector &x_new) const override
        {
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
         bool restrict(const Vector &x, Vector &x_new) const override
        {
            assert(_I);
            x_new = transpose(*_I) * x;

            if(restrict_factor_ != 1.0) {
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
         bool boolean_restrict_or(const Vector &x, Vector &x_new) override
        {
        	assert(false && "implement me");
            // static const Scalar off_diag_tol = std::numeric_limits<Scalar>::epsilon() * 1e6;

            // // x_new = local_zeros(local_size(*_R).get(0));

            // Matrix R_boolean = *_R;
            // R_boolean *= 0.;

            // {
            //     Write<Matrix> w_(R_boolean);

            //     each_read(*_R, [&R_boolean](const SizeType i, const SizeType j, const Scalar value) {
            //         if(std::abs(value) > off_diag_tol) {
            //             R_boolean.set(i, j, 1.);
            //         }
            //     });
            // }

            // x_new = R_boolean * x;

            // ReadAndWrite<Vector> rw_(x_new);

            // auto r = range(x_new);
            // for(auto i = r.begin(); i < r.end(); ++i) {
            //     if(x_new.get(i) > 1.) {
            //         x_new.set(i, 1.);
            //     }
            // }

            // //THIS works for serial (use it once we have a is_parallel and is_serial query)
            // // each_read(*_R, [&x, &x_new](const SizeType i, const SizeType j, const Scalar value) {
            // //     if(x.get(j) != 0. && std::abs(value) > off_diag_tol) {
            // //         x_new.set(i, 1.);
            // //     }
            // // });

            return true;
        }

        void handle_equality_constraints(const Vector &is_constrained) override {
            assert(_I);
            assert(!empty(*_I));
            assert(!empty(is_constrained));
            assert( size(is_constrained).get(0) == size(*_I).get(0) );

            if(empty(is_constrained) || empty(*_I)) return;

            set_zero_rows(*_I, is_constrained);
        }

        /**
         * @brief      Restriction of matrix.
         *
         *             \f$  M_{new} = I^{T} * M  * I  \f$
         * @param[in]  M
         * @param      M_new
         *
         */
        bool restrict(const Matrix &M, Matrix &M_new) const override
        {
            assert(_I);
            M_new =  utopia::ptap(M, *_I);
            return true;
        }


        /**
         * @brief      Projection of vector
         *            \f$  x_{new} = P * x  \f$
         * @param[in]  x
         * @param      x_new
         *
         */
        bool project_down(const Vector &x, Vector &x_new) const override
        {
            assert(_Pr || _I);

            if(!_Pr) {
            	return restrict(x, x_new);
            } else {
            	x_new = *_Pr * x;
        	}
            return true;
        }


        Scalar interpolation_inf_norm() const override
        {
            return norm_infty(*_I);
        }

        Scalar projection_inf_norm() const override
        {
        	m_utopia_warning_once("projection_inf_norm not implemented properly");
            return norm_infty(transpose(*_I));
        }

        Scalar restriction_inf_norm() const override
        {
            return norm_infty(*_Pr);
        }

        protected:
            std::shared_ptr<Matrix> _I;
            std::shared_ptr<Matrix> _Pr;
            Scalar restrict_factor_;
    };

}

#endif //UTOPIA_IP_TRANSFER_HPP

