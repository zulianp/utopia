#ifndef UTOPIA_NONLINEAR_SMOOTHER_HPP
#define UTOPIA_NONLINEAR_SMOOTHER_HPP

#include "utopia_Core.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_Input.hpp"



namespace utopia {

    /**
     * @brief      Nonlinear smoother class
     *  
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector>
    class NonLinearSmoother : virtual public Configurable
    {
        typedef UTOPIA_SCALAR(Vector)                   Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                SizeType;

        public:
        NonLinearSmoother(): _sweeps(1), _relaxation_parameter(1.0)
        { 
            
        }


        virtual void read(Input &in) override
        {
            in.get("relaxation_parameter", _relaxation_parameter);
            in.get("sweeps", _sweeps);

        }


        virtual void print_usage(std::ostream &os) const override
        {
            this->print_param_usage(os, "relaxation_parameter", "double", "Relaxation parameter.", "1"); 
            this->print_param_usage(os, "sweeps", "int", "Number of smoothing steps.", "1"); 
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


        private:
            SizeType     _sweeps;  
            Scalar       _relaxation_parameter; 
};

}

#endif //UTOPIA_NONLINEAR_SMOOTHER_HPP

