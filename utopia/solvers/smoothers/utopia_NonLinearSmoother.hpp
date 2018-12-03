#ifndef UTOPIA_NONLINEAR_SMOOTHER_HPP
#define UTOPIA_NONLINEAR_SMOOTHER_HPP

#include "utopia_Core.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Function.hpp"



namespace utopia {

    /**
     * @brief      Nonlinear smoother class
     *  
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector>
    class NonLinearSmoother 
    {
        typedef UTOPIA_SCALAR(Vector)                   Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                SizeType;

        public:
        NonLinearSmoother(const Parameters params = Parameters()) 
        { 
            set_parameters(params); 
        }


        virtual void set_parameters(const Parameters params) 
        {
            _sweeps = params.pre_smoothing_steps();            
            _relaxation_parameter = params.omega();     
        }


        virtual bool smooth(Function<Matrix, Vector> &fun,  Vector &x, const Vector &rhs) = 0; 


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
         * @brief      Set omega
         *
         * @return     omega  The relaxation parameter.
         */
        virtual bool relaxation_parameter(const Scalar & alpha)
        {
             _relaxation_parameter = alpha; 
             return true; 
        }


        /**
         * @brief      verbose
         *
         */
        virtual void verbose(const bool & verbose)
        {
            _verbose = verbose;
        }

        /**
         * @brief      Verbose
         *
         */
        virtual bool verbose()
        {
            return _verbose; 
        }


        private:
            SizeType     _sweeps;  
            Scalar       _relaxation_parameter; 
            bool         _verbose; 
};

}

#endif //UTOPIA_NONLINEAR_SMOOTHER_HPP

