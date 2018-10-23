#ifndef UTOPIA_SOLVER_QUASI_NEWTON_HPP
#define UTOPIA_SOLVER_QUASI_NEWTON_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_LS_Strategy.hpp"
#include "utopia_HessianApproximations.hpp"
#include "utopia_MatrixFreeSolverInterface.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{

    template<class Matrix, class Vector>
    class QuasiNewton : public NonLinearSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                               Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                            SizeType;
        
        typedef utopia::LSStrategy<Matrix, Vector>                  LSStrategy;
        typedef utopia::HessianApproximation<Matrix, Vector>        HessianApproximation;
        
        typedef utopia::MatrixFreeLinearSolver<Vector>              Solver;
        
    public:

        QuasiNewton(const std::shared_ptr <HessianApproximation> &hessian_approx,
                    const std::shared_ptr <Solver> &linear_solver,
                    const Parameters params = Parameters()):
                    NonLinearSolver<Matrix, Vector>(), 
                    alpha_(1.0),
                    hessian_approx_strategy_(hessian_approx), 
                    mf_linear_solver_(linear_solver)
        {
            set_parameters(params);
        }
        
        bool solve(Function<Matrix, Vector> &fun, Vector &x) override
        {
            using namespace utopia;
            
            Vector g, s, y; 
            
            Scalar g_norm, g0_norm, r_norm=1, s_norm=1;
            SizeType it = 0;
            
            bool converged = false;
            
            fun.gradient(x, g);
            g0_norm = norm2(g);
            g_norm = g0_norm;
            
            if(this->verbose_) {
                this->init_solver("QUASI NEWTON", {" it. ", "|| g ||", "r_norm", "|| p_k || ", "alpha"});
                PrintInfo::print_iter_status(it, {g_norm, r_norm, s_norm});
            }
            it++; 
            
            hessian_approx_strategy_->initialize(fun, x);

            while(!converged)
            {
                if(QuasiLinearSolver<Vector> * lin_solver = dynamic_cast<QuasiLinearSolver<Vector>*>(mf_linear_solver_.get()))
                {
                    auto inverse_action = FunctionOperator<Vector>(hessian_approx_strategy_->get_apply_Hinv()); 
                    lin_solver->solve(inverse_action, -1.0*g, s); 
                }
                else
                {
                    auto multiplication_action = FunctionOperator<Vector>(hessian_approx_strategy_->get_apply_H()); 
                    lin_solver->solve(multiplication_action, -1.0*g, s); 
                }

                if(ls_strategy_) 
                    ls_strategy_->get_alpha(fun, g, x, s, alpha_);     

                s *= alpha_; 
                x+=s; 
                
                y = g; 
                fun.gradient(x, g);

                // norms needed for convergence check
                g_norm = norm2(g);
                r_norm = g_norm/g0_norm;
                s_norm = norm2(s);

                // diff between fresh and old grad...
                y = g - y; 
                hessian_approx_strategy_->update(s, y);

                // print iteration status on every iteration
                if(this->verbose_)
                    PrintInfo::print_iter_status(it, {g_norm, r_norm, s_norm, alpha_});
                
                // check convergence and print interation info
                converged = this->check_convergence(it, g_norm, r_norm, s_norm);
                
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
         * @brief      Sets strategy for computing step-size.
         *
         * @param[in]  strategy  The line-search strategy.
         *
         * @return
         */
        virtual bool set_hessian_approximation_strategy(const std::shared_ptr<HessianApproximation> &strategy)
        {
            hessian_approx_strategy_      = strategy;
            return true;
        }
        
        
    protected:
        Scalar alpha_;                                          /*!< Dumping parameter. */
        std::shared_ptr<LSStrategy> ls_strategy_;               /*!< Strategy used in order to obtain step \f$ \alpha_k \f$ */
        
        std::shared_ptr<HessianApproximation> hessian_approx_strategy_;
        std::shared_ptr<Solver> mf_linear_solver_;                      /*!< Linear solver parameters. */  
        
    };
    
}
#endif //UTOPIA_SOLVER_QUASI_NEWTON_HPP
