#ifndef UTOPIA_MATRIX_TRUNCATED_TRANSFER_HPP
#define UTOPIA_MATRIX_TRUNCATED_TRANSFER_HPP

#include "utopia_Transfer.hpp"

#include <cassert>
#include <cmath>
#include <memory>



     namespace utopia
     {
        /**
         * @brief      The class for transfer operators.
         *
         * @tparam     Matrix
         * @tparam     Vector
         */
        template<class Matrix, class Vector>
        class MatrixTruncatedTransfer final : public Transfer<Matrix, Vector> {
        public:
            using SizeType = typename utopia::Traits<Vector>::SizeType;
            using Scalar = typename utopia::Traits<Vector>::Scalar;

            MatrixTruncatedTransfer(const std::shared_ptr<Matrix> & I)
            {
                assert(I);

                _I = I;
                _R = std::make_shared<Matrix>(transpose(*I));
                // _Pr = _R;

                _R_truncated = std::make_shared<Matrix>(*_R); 
                _I_truncated = std::make_shared<Matrix>(*_I); 
            }

            MatrixTruncatedTransfer(const std::shared_ptr<Matrix> &I, const std::shared_ptr<Matrix> &P):
                    _I(I),
                    _R(std::make_shared<Matrix>(transpose(*I))),
                    _Pr(P)
            {
                assert(I);
                assert(P);

                _R_truncated = std::make_shared<Matrix>(*_R); 
                _I_truncated = std::make_shared<Matrix>(*_I);                 

                // std::cout<<"proper transfer down ... \n";
            }


            MatrixTruncatedTransfer(const std::shared_ptr<Matrix> &I, const std::shared_ptr<Matrix> &R, const std::shared_ptr<Matrix> &P):
                    _I(I),
                    _R(R),
                    _Pr(P)
            {
                assert(I);
                assert(R);
                assert(P);
            }

            ~MatrixTruncatedTransfer() override = default;

            /*=====================================================
                                initialization
            =====================================================*/

            /**
             * @brief      Initialization of interpolation operator.
             *
             * @param[in]  I_in  The interpolation.
             *
             */
            bool I_init(const std::shared_ptr<Matrix> &I_in)
            {
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
            bool R_init(const std::shared_ptr<Matrix> &R_in)
            {
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
            bool IR_init( const std::shared_ptr<Matrix> &I_in, const std::shared_ptr<Matrix> &R_in)
            {
                assert(I_in);
                assert(R_in);

                _I = I_in;
                _R = R_in;

                return true;
            }


            bool scale_transfer(const Matrix &scaling_mat)
            {
                assert(_I);
                assert(_R);

                *_I = scaling_mat* (*_I);
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
             bool interpolate(const Vector &x, Vector &x_new) const override
            {
                assert(_I_truncated);
                x_new = *_I_truncated * x;
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

                assert(_R_truncated);
                x_new = *_R_truncated * x;
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
                static const Scalar off_diag_tol = std::numeric_limits<Scalar>::epsilon() * 1e6;
                if(!_R_truncated){
                    *_R_truncated = *_R; 
                }

                Matrix R_boolean = *_R_truncated;

                each_apply(R_boolean, [](const Scalar value) -> Scalar {
                    if(std::abs(value) > off_diag_tol) {
                        return 1.;
                    } else {
                        return 0.;
                    }
                });

                x_new = R_boolean * x;

                ReadAndWrite<Vector> rw_(x_new);

                auto r = range(x_new);
                for(auto i = r.begin(); i < r.end(); ++i) {
                    if(x_new.get(i) > 1.) {
                        x_new.set(i, 1.);
                    }
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
            bool restrict(const Matrix &M, Matrix &M_new) const override
            {
                assert(_I_truncated);
                M_new =  utopia::ptap(M, *_I_truncated);
                return true;
            }

            /**
             * @brief      Initialization of projection down operator.
             *
             * @param[in]  P_in  The projection operator.
             *
             */
            bool P_init(const std::shared_ptr<Matrix> &P_in)
            {
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
            bool project_down(const Vector &x, Vector &x_new) const override
            {
                assert(_Pr);
                x_new = *_Pr * x;
                return true;
            }

            bool project_down_positive_negative(const Vector &x_pos, const Vector &x_neg, Vector &x_new) override
            {
                if(empty(P_pos_))
                {
                    P_pos_ = *_Pr;
                    chop_smaller_than(P_pos_, 1e-13); 
                }

                if(empty(P_neg_))
                {
                    P_neg_ = (*_Pr); 
                    chop_greater_than(P_neg_, -1e-13); 
                }
                    
                x_new = (P_pos_*x_pos) + (P_neg_*x_neg); 
                return true; 
            }

            const Matrix &I() const
            {
                return *_I;
            }

            const Matrix &R() const
            {
                return *_R;
            }

            Scalar interpolation_inf_norm() const override
            {
                return norm_infty(*_I);
            }

            Scalar projection_inf_norm() const override
            {
                return norm_infty(*_R);
            }

            Scalar restriction_inf_norm() const override
            {
                return norm_infty(*_Pr);
            }

            void truncate_interpolation(const Vector & _eq_active_flg) 
            {
                // to speed up, we should check if constraint was changed between iterations 
                std::vector<SizeType> indices_eq_constraints_; 
                {
                    Read<Vector> r(_eq_active_flg);

                    Range range_w = range(_eq_active_flg);
                    for (SizeType i = range_w.begin(); i != range_w.end(); i++)
                    {
                        if(_eq_active_flg.get(i) == 1.0){
                            indices_eq_constraints_.push_back(i);
                        }
                    }
                }                    

                *_I_truncated = *_I; 
                set_zero_rows(*_I_truncated, indices_eq_constraints_, 0.0);
                *_R_truncated = transpose(*_I_truncated); 
            }

        private:
            std::shared_ptr<Matrix> _I, _R; 

            std::shared_ptr<Matrix> _Pr; // used only for nonlinear mutlilevel solvers
            Matrix P_pos_; // used only for nonlinear mutlilevel solvers
            Matrix P_neg_; // used only for nonlinear mutlilevel solvers

            std::shared_ptr<Matrix> _I_truncated; 
            std::shared_ptr<Matrix> _R_truncated; 
    };

}

#endif //UTOPIA_MATRIX_TRUNCATED_TRANSFER_HPP

