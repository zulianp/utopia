/*
* @Author: Alena Kopanicakova
* @Date:   2017-01-31
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-01-31
*/
#ifndef UTOPIA_SOLVER_INEXACT_NEWTON_HPP
#define UTOPIA_SOLVER_INEXACT_NEWTON_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_LS_Strategy.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{
    /**
     * @brief      The Inexact Newton solver.
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<class Matrix, class Vector>
    class InexactNewton : public NonLinearSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::LSStrategy<Matrix, Vector> LSStrategy; 
        typedef utopia::IterativeSolver<Matrix, Vector> Solver; 
        

    public:
       InexactNewton(   const std::shared_ptr <Solver> &linear_solver = std::shared_ptr<Solver>(),
                        const Parameters params                       = Parameters() ):
                        NonLinearSolver<Matrix, Vector>(linear_solver, params), alpha_(1), hessian_refresh_(4), linear_solver_it_(1)
                        {
                            set_parameters(params);
                            linear_solver->max_it(linear_solver_it_); 
                        }

        bool solve(Function<Matrix, Vector> &fun, Vector &x) override
        {
           using namespace utopia;

            Vector grad, step;
            Matrix hessian;

            Scalar g_norm, g0_norm, r_norm, s_norm;
            SizeType it = 0;

            bool converged = false;

            fun.gradient(x, grad);
            g0_norm = norm2(grad);
            g_norm = g0_norm;

            this->init_solver("INEXACT NEWTON", {" it. ", "|| g ||", "r_norm", "|| p_k || "});

            if(this->verbose_)
                PrintInfo::print_iter_status(it, {g_norm, 1, 0});
        
            while(!converged)
            {
                if(it % hessian_refresh_ == 0)
                    fun.hessian(x, hessian);

                //find direction step
                step = local_zeros(local_size(x));
                this->linear_solve(hessian, grad, step);

                if(ls_strategy_) {
                    
                    ls_strategy_->get_alpha(fun, grad, x, -step, alpha_); 
                    x -= alpha_ * step;

                } else { 
                    //update x
                    if (fabs(alpha_ - 1) < std::numeric_limits<Scalar>::epsilon())
                    {
                        x -= step;
                    }
                    else
                    {
                        x -= alpha_ * step;
                    }
                }

                fun.gradient(x, grad);

                // norms needed for convergence check
                g_norm = norm2(grad);
                r_norm = g_norm/g0_norm;
                s_norm = norm2(step);

                it++;
                // // print iteration status on every iteration
                if(this->verbose_)
                    PrintInfo::print_iter_status(it, {g_norm, r_norm, s_norm});

                // // check convergence and print interation info
                converged = this->check_convergence(it, g_norm, r_norm, s_norm);
            }

            return true;
        }

    virtual void set_parameters(const Parameters params) override
    {
        NonLinearSolver<Matrix, Vector>::set_parameters(params);
        alpha_ = params.alpha();

    }


    /**
     * @brief      Sets strategy for computing step-size. 
     *
     * @param[in]  strategy  The line-search strategy.
     *
     * @return     
     */
    virtual bool set_line_search_strategy(const std::shared_ptr<LSStrategy> &strategy)
    {
      ls_strategy_ = strategy; 
      ls_strategy_->set_parameters(this->parameters());
      return true; 
    }


    /**
     * @return     number of it. to keep same hessian
     */
    virtual const SizeType hessian_refresh()
    {
        return hessian_refresh_; 
    }

    private:
        Scalar alpha_;                                          /*!< Dumping parameter. */
        const SizeType hessian_refresh_;                        /*!< How often refresh hessian. */
        const SizeType linear_solver_it_;                       /*!< How many linear solver iterations required. */
        std::shared_ptr<LSStrategy> ls_strategy_;               /*!< Strategy used in order to obtain step \f$ \alpha_k \f$ */  

    };

}
#endif //UTOPIA_SOLVER_INEXACT_NEWTON_HPP
