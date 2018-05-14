/*
* @Author: alenakopanicakova
* @Date:   2016-04-02
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-07-03
*/

#ifndef UTOPIA_ML_TRANSFER_HPP
#define UTOPIA_ML_TRANSFER_HPP

#include "utopia_Smoother.hpp"
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

        Transfer(const std::shared_ptr<Matrix> & I)//:
                                    // _I(I),
                                    // _R(transpose(I))
        {
            assert(I);

            _I = I; 
            _R = std::make_shared<Matrix>(transpose(*I));
            _Pr = _R;
        }


        // TODO:: missing setter for restriction ... 
        Transfer(const std::shared_ptr<Matrix> &I, const std::shared_ptr<Matrix> &P):
                _I(I),
                _R(std::make_shared<Matrix>(transpose(*I))),
                _Pr(P)
        {
            assert(I);
            assert(P);
            
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
        virtual bool I_init(const std::shared_ptr<Matrix> &I_in)
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
        virtual bool R_init(const std::shared_ptr<Matrix> &R_in)
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
        virtual bool IR_init(
            const std::shared_ptr<Matrix> &I_in, 
            const std::shared_ptr<Matrix> &R_in)
        {
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
        virtual bool interpolate(const Vector &x, Vector &x_new)
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
        virtual bool restrict(const Vector &x, Vector &x_new)
        {

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
        virtual bool boolean_restrict_or(const Vector &x, Vector &x_new)
        {
            static const Scalar off_diag_tol = std::numeric_limits<Scalar>::epsilon() * 1e6;

            // x_new = local_zeros(local_size(*_R).get(0));

            Matrix R_boolean = *_R;
            R_boolean *= 0.;

            {
                Write<Matrix> w_(R_boolean);

                each_read(*_R, [&R_boolean](const SizeType i, const SizeType j, const Scalar value) {
                    if(std::abs(value) > off_diag_tol) {
                        R_boolean.set(i, j, 1.);
                    } 
                });
            }

            x_new = R_boolean * x;

            ReadAndWrite<Vector> rw_(x_new);

            auto r = range(x_new);
            for(auto i = r.begin(); i < r.end(); ++i) {
                if(x_new.get(i) > 1.) {
                    x_new.set(i, 1.);
                } 
            }

            //THIS works for serial (use it once we have a is_parallel and is_serial query)
            // each_read(*_R, [&x, &x_new](const SizeType i, const SizeType j, const Scalar value) {
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
        virtual bool restrict(const Matrix &M, Matrix &M_new) const
        {
            assert(_I);
            M_new =  utopia::ptap(M, *_I);  
            return true; 
        }

        /**
         * @brief      Initialization of projection down operator.
         *
         * @param[in]  P_in  The projection operator. 
         *
         */
        virtual bool P_init(const std::shared_ptr<Matrix> &P_in)
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
        virtual bool project_down(const Vector &x, Vector &x_new)
        {
            assert(_Pr);
            x_new = *_Pr * x;
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

        protected:        
            std::shared_ptr<Matrix> _I, _R; // _P;  
            std::shared_ptr<Matrix> _Pr;
    };

}

#endif //UTOPIA_ML_TRANSFER_HPP

