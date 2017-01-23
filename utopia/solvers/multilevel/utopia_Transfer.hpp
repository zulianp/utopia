/*
* @Author: alenakopanicakova
* @Date:   2016-04-02
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-11-15
*/

#ifndef UTOPIA_ML_TRANSFER_HPP
#define UTOPIA_ML_TRANSFER_HPP
#include "utopia_Smoother.hpp"

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

        Transfer()
        {

        }


        Transfer(const Matrix & I):
                                    _I(I),
                                    _R(transpose(I))
        {

        }

        virtual ~Transfer(){} 
        
        /*=====================================================
                            initialization
        =====================================================*/
        
        /**
         * @brief      Initialization of interpolation operator.
         *
         * @param[in]  I_in  The interpolation. 
         *
         */
        virtual bool I_init(const Matrix &I_in)
        {
            _I = I_in; 
            return true; 
        }

        /**
         * @brief      Initialization of restriction operator.
         *
         * @param[in]  R_in  The restriction. 
         *
         */
        virtual bool R_init(const Matrix &R_in)
        {
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
        virtual bool IR_init(const Matrix &I_in, const Matrix &R_in)
        {
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
        virtual bool interpolate(const Vector &x, Vector &x_new)
        {
            x_new = _I * x; 
            return true; 
        }
        
        /**
         * @brief      Interpolation of matrix. 
         *             \f$  M_{new} = I * M  * I^{T}  \f$
         *
         * @param[in]  M     
         * @param      M_new 
         *
         */
        virtual bool interpolate(const Matrix &M, Matrix &M_new)
        {            
            M_new =  mat_PtAP_product(M, _R); 
            return true; 
        }


        /**
         * @brief      Restriction of vector. 
         *            \f$  x_{new} = R * x  \f$
         * @param[in]  x     
         * @param      x_new 
         *
         */
        virtual bool restrict(const Vector &x, Vector &x_new)
        {
            x_new = _R * x; 
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
        virtual bool restrict(const Matrix &M, Matrix &M_new)
        {
            M_new =  mat_PtAP_product(M, _I);  // petsc implementation of this product is way faster ... 
            return true; 
        }

        /**
         * @brief      Zeros rows of the interpolation operator and saves adjusted I & R on given level. 
         *             Should be called just on the finest level. 
         *
         * @param[in]  I          Interpolation operator.
         * @param[in]  zero_rows  Rows to be zero-ed out. 
         *
         */
        virtual bool apply_truncated_basis_to_interpolation(const std::vector<SizeType>& active_set)
        {
            set_zero_rows(_I, active_set);
            _R = transpose(_I); 
            return true; 
        }


        /**
         * @brief      Return interpolation operator. 
         *
         */
        Matrix I()
        {
            return _I; 
        }


    private:        
        Matrix _I, _R;  


    };

}

#endif //UTOPIA_ML_TRANSFER_HPP

