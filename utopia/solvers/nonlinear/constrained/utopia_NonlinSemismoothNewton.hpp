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
#include "utopia_BoxConstraints.hpp"  
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
        typedef utopia::BoxConstraints<Vector>      BoxConstraints;

    public:
       NonlinSemismoothNewton(  const std::shared_ptr <Solver> &linear_solver   = std::shared_ptr<Solver>(),
                                const Parameters   params                       = Parameters())
                                                    : NonLinearSolver<Matrix, Vector>(linear_solver, params)
        {  

        }

        bool solve(Function<Matrix, Vector> & fun, Vector & x_new) override 
        {
             using namespace utopia;
    
            Scalar c = 1;
            SizeType iterations = 1;

            const SizeType N = x_new.size().get(0);
            
            Vector lambda = zeros(N);
            Vector Ginvg, d, g;
            Vector x_old = x_new;
            
            Matrix Ac, Ic, M;
            Matrix G = identity(N,N);

            Matrix Hessian;
            Vector grad;

            Vector upbo; 
            if(constraints_->has_upper_bound())
                upbo = *constraints_->upper_bound(); 
            else
            {
                // TODO:: fix 
                std::cout<<"NonlinSemismoothNewton does not support other types at the moment.... \n"; 
            }

            this->linear_solve(G, upbo, Ginvg);

            bool converged = false;
            this->init_solver("NON-LINEAR - SEMISMOOTH NEWTON METHOD", {" it. ", " err "});

            while(!converged)
            {
                // this is super expensive 
                d = lambda + c * (G * x_new - upbo); 

                //!
                //! active, inactive constraints 
                Ac = zeros(N, N);
                Ic = zeros(N, N);
                {
                    Read<Vector> r (d);
                    Write<Matrix> w_Ac (Ac);
                    Write<Matrix> w_Ic (Ic);
                    Range rr = rowRange(Ac);    
                    for ( SizeType i = rr.begin(); i != rr.end(); i ++){
                        if(d.get(i) > 0){    
                                Ac.set(i,i, 1.0);
                        }
                       else {
                                Ic.set(i,i, 1.0);
                        }
                    }
                }
                
                // if(b_range.begin() == 0) b.set(0, 0);
                // if(b_range.end() == N) b.set(N - 1, 0);


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


        virtual bool set_box_constraints(const std::shared_ptr<BoxConstraints> & box)
        {
          constraints_ = box; 
          return true; 
        }

        virtual std::shared_ptr<BoxConstraints> get_box_constraints() const
        {
          return constraints_; 
        }


        virtual void set_parameters(const Parameters params) override
        {
            NonLinearSolver<Matrix, Vector>::set_parameters(params);
        }
    
    private:
        std::shared_ptr<BoxConstraints> constraints_; 
    };

}
#endif //UTOPIA_SOLVER_NONLINSEMISMOOTH_NEWTON_HPP
