/*
* @Author: alenakopanicakova
* @Date:   2016-04-02
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-07-03
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

        Transfer(const std::shared_ptr <Matrix> & I)//:
                                    // _I(I),
                                    // _R(transpose(I))
        {
            _I = I; 
            _R = std::make_shared<Matrix>(transpose(*I)); 
            _P = _R; 
        }

        

        Transfer(const std::shared_ptr <Matrix> & I, const std::shared_ptr <Matrix> & P):
                _I(I),
                _R(std::make_shared<Matrix>(transpose(*I))),
                _P(P)
        {
         std::cout<<"proper transfer down ... \n"; 
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
        virtual bool I_init(const std::shared_ptr <Matrix> &I_in)
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
        virtual bool R_init(const std::shared_ptr <Matrix> &R_in)
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
        virtual bool IR_init(const std::shared_ptr <Matrix> &I_in, const std::shared_ptr <Matrix> &R_in)
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
            x_new = *_I * x; 
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
        // virtual bool interpolate(const std::shared_ptr < Matrix> &M, std::shared_ptr <Matrix> &M_new)
        // {            
        //     M_new =  std::make_shared<Matrix>(mat_PtAP_product(*M, *_R)); 
        //     return true; 
        // }

        /**
         * @brief      Interpolation of matrix. 
         *             \f$  M_{new} = I * M  * I^{T}  \f$
         *
         * @param[in]  M     
         * @param      M_new 
         *
         */
        // virtual bool interpolate(const std::shared_ptr <const Matrix> &M, std::shared_ptr <Matrix> &M_new)
        // {            
        //     M_new =  std::make_shared<Matrix>(mat_PtAP_product(*M, *_R)); 
        //     return true; 
        // }

        /**
         * @brief      Restriction of vector. 
         *            \f$  x_{new} = R * x  \f$
         * @param[in]  x     
         * @param      x_new 
         *
         */
        virtual bool restrict(const Vector &x, Vector &x_new)
        {
            x_new = *_R * x; 
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
        // virtual bool restrict(const std::shared_ptr < Matrix> &M,  std::shared_ptr <Matrix> &M_new)
        // {
        //     M_new =  std::make_shared<Matrix>(mat_PtAP_product(*M, *_I));  
        //     return true; 
        // }

        /**
         * @brief      Restriction of matrix.
         *             
         *             \f$  M_{new} = I^{T} * M  * I  \f$
         * @param[in]  M     
         * @param      M_new 
         *
         */
        // virtual bool restrict(const std::shared_ptr <const  Matrix> &M,  std::shared_ptr <  Matrix> &M_new)
        // {
        //     M_new =  std::make_shared<Matrix>(mat_PtAP_product(*M, *_I));  
        //     return true; 
        // }


        /**
         * @brief      Restriction of matrix.
         *             
         *             \f$  M_{new} = I^{T} * M  * I  \f$
         * @param[in]  M     
         * @param      M_new 
         *
         */
        virtual bool restrict(const Matrix &M, Matrix &M_new) const
        {
            M_new =  mat_PtAP_product(M, *_I);  
            return true; 
        }



        /**
         * @brief      Initialization of projection down operator.
         *
         * @param[in]  P_in  The projection operator. 
         *
         */
        virtual bool P_init(const std::shared_ptr <Matrix> &P_in)
        {
            _P = P_in; 
            return true; 
        }



        /**
         * @brief      Projection of vector 
         *            \f$  x_{new} = P * x  \f$
         * @param[in]  x     
         * @param      x_new 
         *
         */
        virtual bool project_down(const Vector &x, Vector &x_new)
        {
            x_new = *_P * x; 
            return true; 
        }


        protected:        
            std::shared_ptr <Matrix> _I, _R; // _P;  
            std::shared_ptr <Matrix>  _P;  


    };

}

#endif //UTOPIA_ML_TRANSFER_HPP

