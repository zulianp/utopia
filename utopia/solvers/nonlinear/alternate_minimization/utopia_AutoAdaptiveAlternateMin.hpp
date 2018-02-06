/*
* @Author: kopanicakova
* @Date:   2018-02-02 10:21:12
* @Last Modified by:   kopanicakova
* @Last Modified time: 2018-02-06
*/
#ifndef UTOPIA_AUTO_ADAPTIVE_ALTERNATE_MINIMIZATION_HPP
#define UTOPIA_AUTO_ADAPTIVE_ALTERNATE_MINIMIZATION_HPP

#include "utopia_Core.hpp" 
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"

#include <iomanip>
#include <limits>
#include <cmath>


namespace utopia
{
    template<class Matrix, class Vector, class FunctionType>
    class AutoAdaptiveAlternateMin : public AlternateMinimization<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::NonLinearSolver<Matrix, Vector>     NonlinearSolver;

    public:
       AutoAdaptiveAlternateMin(   const std::shared_ptr<NonlinearSolver> &nl_solver_master = std::shared_ptr<NonlinearSolver>(),
                                    const std::shared_ptr<NonlinearSolver> &nl_solver_slave = std::shared_ptr<NonlinearSolver>(), 
                                    const Parameters params                                 = Parameters()):
                                    AlternateMinimization<Matrix, Vector>(nl_solver_master, nl_solver_slave, params)
        {
            set_parameters(params); 
        }


        bool solve( AlternateMinLevel<Matrix, Vector,FunctionType > level, 
                    Vector &x_mono, 
                    Vector &x_stag1, 
                    Vector &x_stag2) 
        {
            using namespace utopia;

            // // ------------------------ monolithic  --------------------
            std::cout<<"------------- monolithic initial ---------------------- \n"; 
            Vector x_mono_test = x_mono; 
            this->_nl_solver_master->solve(level.fun_monolithic(), x_mono_test); 




            // ------------------------ staggered --------------------
            std::cout<<"------------- staggered ---------------------- \n"; 
            Scalar energy_0, energy_prev, energy_cur; 

            level.fun_stag1().value(x_stag1, energy_0); 
            energy_prev = energy_0; 


            Scalar beta    = 9e9; 
            SizeType it     = 1; 

            while(beta > this->get_energy_slope_tol() && it < this->get_num_alternate_steps() )
            {
                this->_nl_solver_master->max_it(1); 
                this->_nl_solver_master->solve(level.fun_stag1(), x_stag1); 
                level.transfer_from_stag1_to_stag2(x_stag1); 

                this->_nl_solver_slave->solve(level.fun_stag2(), x_stag2); 
                level.transfer_from_stag2_to_stag1(x_stag2); 


                level.fun_stag1().value(x_stag1, energy_cur); 

                Scalar val = ((energy_prev - energy_cur)/ (energy_0 - energy_cur)) * it; 
                beta = std::atan (val) * 180.0 / M_PI;

                energy_prev = energy_cur; 
                it++; 


                level.transfer_from_stag1_to_monolithic(x_stag1, x_mono); 
                level.transfer_from_stag2_to_monolithic(x_stag2, x_mono); 

                // ------------------------ monolithic  --------------------
                std::cout<<"------------- monolithic ---------------------- \n"; 
                x_mono_test = x_mono; 
                this->_nl_solver_master->max_it(50); 
                this->_nl_solver_master->solve(level.fun_monolithic(), x_mono_test); 

                // from monolithic to staggered ... 

            }


            level.transfer_from_stag1_to_monolithic(x_stag1, x_mono); 
            level.transfer_from_stag2_to_monolithic(x_stag2, x_mono); 

            // // ------------------------ monolithic  --------------------
            // std::cout<<"------------- monolithic ---------------------- \n"; 
            // x_mono_test = x_mono; 
            // this->_nl_solver_master->solve(level.monolithic, x_mono_test); 


            // Scalar diff_mono = norm2(x_mono_test - x_mono); 
            // static Scalar time_step = 1; 

            // std::cout<<"------ time_step: " << time_step << "      diff_mono:   "<< diff_mono <<  "  \n"; 
            // time_step++; 


            // stopping criterium: based on norm of global-gradient ??? 


            return true;
        }


        virtual void set_parameters(const Parameters params)
        {
            AlternateMinimization<Matrix, Vector>::set_parameters(params);
        }





    };

}
#endif //UTOPIA_AUTO_ADAPTIVE_ALTERNATE_MINIMIZATION_HPP
