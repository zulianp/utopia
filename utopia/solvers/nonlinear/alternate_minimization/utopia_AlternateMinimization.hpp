/*
* @Author: Alena Kopanicakova
* @Date:    2018-01-25
* @Last Modified by:   kopanicakova
* @Last Modified time: 2018-01-28
*/
#ifndef UTOPIA_SOLVER_ALTERNATE_MINIMIZATION_HPP
#define UTOPIA_SOLVER_ALTERNATE_MINIMIZATION_HPP

#include "utopia_Core.hpp" 
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"

#include <iomanip>
#include <limits>
#include <cmath>


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
                    std::function<void(const Vector &)> transfer_from_master, 
                    std::function<void(const Vector &)> transfer_to_master) 
        {
            using namespace utopia;

            Scalar energy_1, energy_2, energy_0, energy_prev, energy_cur; 

            fun_master.value(x_master, energy_1); 
            fun_slave.value(x_slave, energy_2); 

            energy_0 = energy_1 + energy_2; 
            energy_prev = energy_0; 

            CSVWriter writer; 
            writer.open_file("/Users/alenakopanicakova/Desktop/tex_files/papers/multilevel_for_PF/results/alternate_2D_tension_beta_test.csv"); 
            // writer.open_file("alternate_2D_tension_beta_test.csv"); 
            // writer.write_table_row<std::string>({"it", "beta"}); 

            Scalar beta    = 9e9; 
            SizeType it     = 1; 

            while(beta > _energy_slope_tol && it < _num_alternate_steps)
            {
               _nl_solver_master->solve(fun_master, x_master); 
               transfer_from_master(x_master); 

               _nl_solver_slave->solve(fun_slave, x_slave); 
               transfer_to_master(x_slave); 

                fun_master.value(x_master, energy_1); 
                fun_slave.value(x_slave, energy_2); 
                energy_cur = energy_1 + energy_2; 

                Scalar val = ((energy_prev - energy_cur)/ (energy_0 - energy_cur)) * it; 
                beta = std::atan (val) * 180.0 / M_PI;
    
                // writer.write_table_row<Scalar>({Scalar(it), beta}); 

                energy_prev = energy_cur; 
                it++; 
            }

            writer.write_table_row<SizeType>({(it-1)}); 
            writer.close_file(); 
            return true;
        }


        virtual void set_parameters(const Parameters params)
        {
            _num_alternate_steps    = params.num_alternate_steps();
            _energy_slope_tol       = params.energy_slope_tol(); 

        }

        virtual void set_num_alternate_steps(const SizeType & num_alternate_steps_in)
        {
            _num_alternate_steps = num_alternate_steps_in; 
        }

        virtual SizeType get_num_alternate_steps()
        {
            return _num_alternate_steps;
        }
        

        virtual void set_energy_slope_tol(const Scalar & energy_slope_tol)
        {
            _energy_slope_tol = energy_slope_tol; 
        }

        virtual Scalar get_energy_slope_tol()
        {
            return _energy_slope_tol;
        }


    private:

        SizeType _num_alternate_steps; 
        Scalar   _energy_slope_tol;                 // in degrees 

        std::shared_ptr<NonlinearSolver>             _nl_solver_master;  
        std::shared_ptr<NonlinearSolver>             _nl_solver_slave;  

    };

}
#endif //UTOPIA_SOLVER_ALTERNATE_MINIMIZATION_HPP
