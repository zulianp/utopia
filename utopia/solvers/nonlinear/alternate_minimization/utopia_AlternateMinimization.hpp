#ifndef UTOPIA_SOLVER_ALTERNATE_MINIMIZATION_HPP
#define UTOPIA_SOLVER_ALTERNATE_MINIMIZATION_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{
    template<class Matrix, class Vector>
    class AlternateMinimization 
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::NonLinearSolver<Matrix, Vector>     NonlinearSolver;

    public:
       AlternateMinimization(   const std::shared_ptr<NonlinearSolver> &nl_solver_master = std::shared_ptr<NonlinearSolver>(),
                                const std::shared_ptr<NonlinearSolver> &nl_solver_slave = std::shared_ptr<NonlinearSolver>(), 
                                const Parameters params                                 = Parameters()):
                            _nl_solver_master(nl_solver_master), 
                            _nl_solver_slave(nl_solver_slave) 
        {
            set_parameters(params); 
        }




        bool solve( Function<Matrix, Vector> &fun_master, Function<Matrix, Vector> &fun_slave, 
                    Vector &x_master, Vector &x_slave, 
                    std::function<void(const Vector &)> multiapp_transfer_from_master, 
                    std::function<void(const Vector &)> multiapp_transfer_to_master) 
        {
           using namespace utopia;

           for(auto it=0; it < _num_alternate_steps; it++)
           {
               _nl_solver_master->solve(fun_master, x_master); 
               multiapp_transfer_from_master(x_master); 


               _nl_solver_slave->solve(fun_slave, x_slave); 
               multiapp_transfer_to_master(x_slave); 
           }

            return true;
        }


        virtual void set_parameters(const Parameters params)
        {
            _num_alternate_steps = params.num_alternate_steps();
        }

        virtual void set_num_alternate_steps(const SizeType & num_alternate_steps_in)
        {
            _num_alternate_steps = num_alternate_steps_in; 
        }

        virtual SizeType set_num_alternate_steps()
        {
            return _num_alternate_steps;
        }


    private:

        SizeType _num_alternate_steps; 

        std::shared_ptr<NonlinearSolver>             _nl_solver_master;  
        std::shared_ptr<NonlinearSolver>             _nl_solver_slave;  

    };

}
#endif //UTOPIA_SOLVER_ALTERNATE_MINIMIZATION_HPP
