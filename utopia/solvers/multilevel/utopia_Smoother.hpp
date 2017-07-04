/*
* @Author: alenakopanicakova
* @Date:   2016-04-04
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-07-02
*/

#ifndef UTOPIA_SMOOTHER_HPP
#define UTOPIA_SMOOTHER_HPP
#include "utopia_Core.hpp"
#include "utopia_Parameters.hpp"    
#include <iomanip>
#include "utopia_Function.hpp"



     namespace utopia 
     {
        template<class Matrix, class Vector>
        class Smoother
        {
            typedef UTOPIA_SCALAR(Vector)           Scalar;
            typedef UTOPIA_SIZE_TYPE(Vector)        SizeType;


        public:


        /**
         * @brief      Base class for smoothers. 
         */
        Smoother() 
        {

        }

        virtual ~Smoother() {}

        virtual void set_parameters(const Parameters params)
        {
            _sweeps = params.pre_smoothing_steps();            
            _relaxation_parameter = params.omega();            

        }

        /**
         * @brief      Single sweep. Function needs to be provided by actual smoothers.
         * @return    
         */
        virtual bool smooth(const Matrix &A, const Vector &rhs, Vector &x) = 0; 



        /**
         * @brief      Quick interface for smoothing with projecting constraints.  
         */
        virtual bool nonlinear_smooth(const Matrix &/*A*/, const Vector &/*rhs*/, const Vector& /*ub*/, const Vector& /*lb*/, Vector &/*x*/, std::vector<SizeType>& /*zero_rows*/){ return 0; }

        /**
         * @brief      Get number of sweeps.
         *
         * @return     
         */
        virtual SizeType sweeps()
        {
            return _sweeps; 
        }


        /**
         * @brief      Set the sweeps.
         *
         * @param[in]  sweeps   The number of sweeps. 
         *
         * @return    
         */
        virtual bool sweeps(const SizeType & sweeps_in)
        {
            _sweeps = sweeps_in;
            return true; 
        }


        /**
         * @brief      Set omega.
         *
         * @return     omega  The relaxation parameter.
         */
        virtual Scalar relaxation_parameter()
        {
            return _relaxation_parameter; 
        }


        /**
         * @brief      Set omega.
         *
         * @return     omega  The relaxation parameter.
         */
        virtual bool relaxation_parameter(const Scalar & relaxation_parameter)
        {
             _relaxation_parameter = relaxation_parameter; 
             return true; 
        }

    private:
        SizeType     _sweeps;  
        Scalar       _relaxation_parameter; 
};

}

#endif //UTOPIA_SMOOTHER_HPP

