#ifndef UTOPIA_SOLVER_INEXACT_NEWTON_HPP
#define UTOPIA_SOLVER_INEXACT_NEWTON_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_LS_Strategy.hpp"
#include "utopia_HessianApproximations.hpp"

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
        typedef UTOPIA_SCALAR(Vector)                           Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;
        
        typedef utopia::LSStrategy<Matrix, Vector>              LSStrategy;
        typedef utopia::HessianApproximation<Matrix, Vector>    HessianApproximation;
        typedef utopia::IterativeSolver<Matrix, Vector>         Solver;
        
        
    public:
        InexactNewton(     const std::shared_ptr <Solver> &linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector> >(),
                      const Parameters params                                         = Parameters()
                      ):
        NonLinearSolver<Matrix, Vector>(linear_solver, params),alpha_(1),
        hessian_refresh_(4),
        linear_solver_it_(1),
        _has_hessian_approx_strategy(false)
        {
            set_parameters(params);
            linear_solver->max_it(linear_solver_it_);
        }
        
        bool solve(Function<Matrix, Vector> &fun, Vector &x) override
        {
            using namespace utopia;
            
            Vector grad, step;
            Matrix hessian;
            
            Scalar g_norm, g0_norm, r_norm=1, s_norm=1;
            SizeType it = 0;
            
            bool converged = false;
            
            fun.gradient(x, grad);
            g0_norm = norm2(grad);
            g_norm = g0_norm;
            
            if(this->verbose_) {
                this->init_solver("INEXACT NEWTON", {" it. ", "|| g ||", "r_norm", "|| p_k || "});
            }
            
            if(has_hessian_approx()) {
                hessian_approx_strategy_->initialize(fun, x, hessian);
            } else {
                fun.hessian(x, hessian);
            }
            
            
            while(!converged)
            {
                // print iteration status on every iteration
                if(this->verbose_)
                    PrintInfo::print_iter_status(it, {g_norm, r_norm, s_norm});
                
                // check convergence and print interation info
                converged = this->check_convergence(it, g_norm, r_norm, s_norm);
                
                if(converged) {
                    return true;
                }
                
                //find direction step
                step = local_zeros(local_size(x));
                this->linear_solve(hessian, -1* grad, step);
                
                if(ls_strategy_) {
                    
                    ls_strategy_->get_alpha(fun, grad, x, step, alpha_);
                    x += alpha_ * step;
                    
                } else {
                    //update x
                    if (fabs(alpha_ - 1) < std::numeric_limits<Scalar>::epsilon())
                    {
                        x += step;
                    }
                    else
                    {
                        x += alpha_ * step;
                    }
                }
                
                if(has_hessian_approx())
                    hessian_approx_strategy_->approximate_hessian(fun, x, step, hessian,  grad);
                else
                    fun.hessian(x, hessian);
                
                // norms needed for convergence check
                g_norm = norm2(grad);
                r_norm = g_norm/g0_norm;
                s_norm = norm2(step);
                
                it++;
            }

            this->print_statistics(it); 
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
        
        /**
         * @brief      Sets strategy for computing step-size.
         *
         * @param[in]  strategy  The line-search strategy.
         *
         * @return
         */
        virtual bool set_hessian_approximation_strategy(const std::shared_ptr<HessianApproximation> &strategy)
        {
            hessian_approx_strategy_      = strategy;
            _has_hessian_approx_strategy  = true;
            return true;
        }
        
        
        virtual bool has_hessian_approx()
        {
            return _has_hessian_approx_strategy;
        }
        
        
        
    private:
        Scalar alpha_;                                          /*!< Dumping parameter. */
        const SizeType hessian_refresh_;                        /*!< How often refresh hessian. */
        const SizeType linear_solver_it_;                       /*!< How many linear solver iterations required. */
        std::shared_ptr<LSStrategy> ls_strategy_;               /*!< Strategy used in order to obtain step \f$ \alpha_k \f$ */
        
        std::shared_ptr<HessianApproximation> hessian_approx_strategy_;
        
        bool _has_hessian_approx_strategy;
        
    };
    
}
#endif //UTOPIA_SOLVER_INEXACT_NEWTON_HPP
