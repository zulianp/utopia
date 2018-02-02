/*
* @Author: kopanicakova
* @Date:   2018-02-02 10:21:12
* @Last Modified by:   kopanicakova
* @Last Modified time: 2018-02-02 10:20:19
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

            std::cout<<"--------- new solve test -------------- \n"; 


            this->_nl_solver_master->solve(level.monolithic, x_mono); 

            level.from_monolithic_to_stag1(x_mono, x_stag1); 
            

            return true;
        }


        virtual void set_parameters(const Parameters params)
        {
            AlternateMinimization<Matrix, Vector>::set_parameters(params);
        }




    };

}
#endif //UTOPIA_AUTO_ADAPTIVE_ALTERNATE_MINIMIZATION_HPP
