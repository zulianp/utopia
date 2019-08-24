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

        bool interpolate(const Vector &x, Vector &x_new) const override
        {
            x_new = x;
            return true;
        }


        bool restrict(const Vector &x, Vector &x_new) const override
        {
            x_new = x;
            return true;
        }


        bool boolean_restrict_or(const Vector &x, Vector &x_new) override
        {
            x_new = x;
            return true;
        }


        bool restrict(const Matrix &M, Matrix &M_new) const override
        {
            M_new = M;
            return true;
        }


        bool project_down(const Vector &x, Vector &x_new) const override
        {
            x_new = x;
            return true;
        }

        bool project_down_positive_negative(const Vector &x_pos, const Vector & /*x_neg*/, Vector &x_new) override
        {
            // as all elements of the identity matrix are positive 
            x_new = x_pos; 
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

