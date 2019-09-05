#ifndef UTOPIA_ML_TRANSFER_HPP
#define UTOPIA_ML_TRANSFER_HPP

#include "utopia_Traits.hpp"

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
        class Transfer
        {
            typedef UTOPIA_SCALAR(Vector)    Scalar;
            typedef UTOPIA_SIZE_TYPE(Vector) SizeType;


        public:


        virtual ~Transfer(){}

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
        virtual bool interpolate(const Vector &x, Vector &x_new) const = 0;

        /**
         * @brief      Restriction of vector.
         *            \f$  x_{new} = R * x  \f$
         * @param[in]  x
         * @param      x_new
         *
         */
        virtual bool restrict(const Vector &x, Vector &x_new) const = 0;

        /**
         * @brief      Restriction of vector based on booleans.
         * values have to be exact 1. and 0.
         *            \f$  x_{new} = R * x  \f$
         * @param[in]  x
         * @param      x_new
         *
         */
        virtual bool boolean_restrict_or(const Vector &x, Vector &x_new) = 0;

        /**
         * @brief      Restriction of matrix.
         *
         *             \f$  M_{new} = I^{T} * M  * I  \f$
         * @param[in]  M
         * @param      M_new
         *
         */
        virtual bool restrict(const Matrix &M, Matrix &M_new) const = 0;

        /**
         * @brief      Projection of vector
         *            \f$  x_{new} = P * x  \f$
         * @param[in]  x
         * @param      x_new
         *
         */
        virtual bool project_down(const Vector &x, Vector &x_new) const = 0;

        /**
         * @brief      Splits projection operator into positive and negative part and then transfers them separatelly, so we end up with
         *            \f$  x_{new} = (P_+ * x_+) + (P_- * x_-)  \f$
         *            
         * @param[in]  x
         * @param      x_new
         *
         */
        virtual bool project_down_positive_negative(const Vector &x_pos, const Vector &x_neg, Vector &x_new) = 0;



        //FIXME
        virtual Scalar interpolation_inf_norm() const = 0;
        virtual Scalar projection_inf_norm()    const = 0;
        virtual Scalar restriction_inf_norm()   const = 0;

        virtual void handle_equality_constraints(const Vector &) {}
    };

}

#endif //UTOPIA_ML_TRANSFER_HPP

