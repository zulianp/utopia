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
       AlternateMinimization(   const std::shared_ptr<NonlinearSolver> &nl_solver1 = std::shared_ptr<NonlinearSolver>(),
                                const std::shared_ptr<NonlinearSolver> &nl_solver2 = std::shared_ptr<NonlinearSolver>()):
                            _nl_solver1(nl_solver1), 
                            _nl_solver2(nl_solver2) 
        {
            //set_parameters(params);
        }

        bool solve(Function<Matrix, Vector> &fun, Vector &x) override
        {
           using namespace utopia;

         
           std::cout<<"-------- alternate minimization solve ----------- \n"; 


            return true;
        }

        // virtual void set_parameters(const Parameters params) override
        // {
        //     NonLinearSolver<Matrix, Vector>::set_parameters(params);
        //     alpha_ = params.alpha();

        // }


    private:

        std::shared_ptr<NonlinearSolver>             _nl_solver1;  
        std::shared_ptr<NonlinearSolver>             _nl_solver2;  

    };

}
#endif //UTOPIA_SOLVER_ALTERNATE_MINIMIZATION_HPP
