#ifndef UTOPIA_IDENTITY_TRANSFER_HPP
#define UTOPIA_IDENTITY_TRANSFER_HPP

#include "utopia_Transfer.hpp"
#include "utopia_Temp.hpp"

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
        class IdentityTransfer final : public Transfer<Matrix, Vector>
        {
            typedef UTOPIA_SCALAR(Vector)    Scalar;

        public:

        IdentityTransfer()
        {

        }


         ~IdentityTransfer(){}


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
            x_new = x; 
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
            x_new = x; 
            return true;
        }

        void handle_equality_constraints(const Vector &is_constrained) override 
        {
            // TO BE thought about... but most likely, we do not need to do anything...
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
            M_new =  M; 
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
            x_new = x;
            return true;
        }


        Scalar interpolation_inf_norm() const override
        {
            return 1.0; 
        }

        Scalar projection_inf_norm() const override
        {
            return 1.0; 
        }

        Scalar restriction_inf_norm() const override
        {
            return 1.0; 
        }


    };

}

#endif //UTOPIA_IDENTITY_TRANSFER_HPP

