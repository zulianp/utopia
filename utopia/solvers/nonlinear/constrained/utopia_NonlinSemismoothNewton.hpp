/*! \file utopia_NonlinSemismoothNewton.hpp
    Nonlinear Semismooth Newton method 
    Created by Alena Kopanicakova 
*/

#ifndef UTOPIA_SOLVER_NONLINSEMISMOOTH_NEWTON_HPP
#define UTOPIA_SOLVER_NONLINSEMISMOOTH_NEWTON_HPP

#include "utopia_NonLinearSolver.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_FunctionBoxConstrained.hpp"  
// #include "utopia_BoxConstaints.hpp"  
#include "utopia_Core.hpp"
#include <vector>

namespace utopia 
{   
    /**
     * @brief      Nonlinear Semi-smooth solver. 
     *
     * @tparam     Matrix 
     * @tparam     Vector 
     */
    template<class Matrix, class Vector>
    class NonlinSemismoothNewton : public NonLinearSolver<Matrix, Vector> 
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef typename NonLinearSolver<Matrix, Vector>::Solver Solver;

    public:
       NonlinSemismoothNewton(  const std::shared_ptr <Solver> &linear_solver   = std::shared_ptr<Solver>(),
        const Parameters   params                       = Parameters())
       : NonLinearSolver<Matrix, Vector>(linear_solver, params)
       {  }

        // void set_box_constaints(const BoxConstaints<Vector> &box_constaints)
        // {
        //     box_constaints_ = box_constaints;
        // }

        //! solve Nonlin - Semismooth Newton
        /*!
          \param fun function
          \param x_new - initial guess/result
          \param A - laplacian
          \param B - boundary operator
          \param tol 
          \param maxIter 
        */
        // template<class Function>
       bool solve(const FunctionBoxConstrained<Matrix, Vector> &fun, Vector &x_new) 
       {
         using namespace utopia;

         Scalar c = 1;
         SizeType iterations = 1;

         const SizeType local_N = local_size(x_new).get(0);

         const SizeType N = x_new.size().get(0);

         Vector lambda = local_zeros(local_N);
         Vector Ginvg, d, g;
         Vector x_old = x_new;

         Matrix Ac, Ic, M;
         Matrix G = local_identity(local_N, local_N);

         Matrix Hessian;
         Vector grad;

         Vector upbo; 
         fun.upper_bound(upbo); 

         this->linear_solve(G, upbo, Ginvg);

         bool converged = false;
         this->init_solver("NON-LINEAR - SEMISMOOTH NEWTON METHOD", { " it. ", " err " });

         while(!converged) {
            // this is super expensive 
            d = lambda + c * (G * x_new - upbo); 

            //! active, inactive constraints 
            if(is_sparse<Matrix>::value) {
                Ac = local_sparse(local_N, local_N, 1);
                Ic = local_sparse(local_N, local_N, 1);
            } else {
                Ac = local_zeros({ local_N, local_N });
                Ic = local_zeros({ local_N, local_N });
            }
            
            {
                Read<Vector> r (d);
                Write<Matrix> w_Ac (Ac);
                Write<Matrix> w_Ic (Ic);
                
                const Range rr = row_range(Ac);    

                for ( SizeType i = rr.begin(); i != rr.end(); i ++){
                    if(d.get(i) > 0){    
                        Ac.set(i, i, 1.0);
                    }
                    else {
                        Ic.set(i, i, 1.0);
                    }
                }
            }

            fun.gradient(x_new, grad);
            fun.hessian(x_new, Hessian);

            g = -1 * grad + Hessian * x_new;

            M = Ac + Ic * Hessian;
            this->linear_solve(M, (Ic * g + Ac * Ginvg), x_new);
            lambda = (g - Hessian * x_new);

            Scalar err = norm2(x_new - x_old);

                // print iteration status on every iteration 
            if(this->verbose_)
                PrintInfo::print_iter_status(iterations, {err}); 


                // check convergence and print interation info
            converged = this->check_convergence(iterations, err, 1, 1);
            x_old = x_new;

            iterations++;
        }
        return true;
    }

    bool solve(Function<Matrix, Vector> & , Vector & ) override 
    {
        return true;
    }


private:
        // BoxConstaints<Vector> box_constaints_;
};

}
#endif //UTOPIA_SOLVER_NONLINSEMISMOOTH_NEWTON_HPP
