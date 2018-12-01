#ifndef UTOPIA_QUASI_NEWTON_BOUND_HPP
#define UTOPIA_QUASI_NEWTON_BOUND_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_LS_Strategy.hpp"
#include "utopia_HessianApproximations.hpp"
#include "utopia_VariableBoundSolverInterface.hpp"
#include "utopia_QPSolver.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{
    template<class Vector>
    class QuasiNewtonBound :    public MatrixFreeNonLinearSolver<Vector>, 
                                public VariableBoundSolverInterface<Vector>,
                                public QuasiNewtonBase<Vector>

    {
        typedef UTOPIA_SCALAR(Vector)                           Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;
        
        typedef utopia::LSStrategy<Vector>                      LSStrategy;
        typedef utopia::HessianApproximation<Vector>            HessianApproximation;

        typedef utopia::MatrixFreeLinearSolver<Vector>          LinSolver;
        // typedef utopia::QPSolver<Vector>                      QPSolver;
        
        
    public:
        QuasiNewtonBound(   const std::shared_ptr <HessianApproximation> &hessian_approx,
                            const std::shared_ptr <LinSolver> &linear_solver,
                            const Parameters params = Parameters()):
                            MatrixFreeNonLinearSolver<Vector>(params), 
                            QuasiNewtonBase<Vector>(hessian_approx, linear_solver), 
                            alpha_(1.0)
        {
            set_parameters(params);
        }
        
        bool solve(FunctionBase<Vector> &fun, Vector &x) override
        {
            using namespace utopia;
            
            Vector g, s, y; 
            
            Scalar g_norm, g0_norm, s_norm=1;
            SizeType it = 0;
            
            bool converged = false;

            this->make_iterate_feasible(x); 
            
            fun.gradient(x, g);
            g0_norm = norm2(g);
            g_norm = g0_norm;
            
            if(this->verbose_) {
                this->init_solver("QUASI NEWTON BOUND", {" it. ", "|| g ||", "E", "|| p_k || ", "alpha"});
                PrintInfo::print_iter_status(it, {g_norm, s_norm});
            }
            it++; 
            
            this->initialize_approximation(); 

            while(!converged)
            {
                // if(TRBoxSubproblem * tr_subproblem = dynamic_cast<TRBoxSubproblem*>(this->linear_solver().get()))
                // {
                //     auto box = this->get_box_constraints();
                //     auto multiplication_action = FunctionOperator<Vector>(this->get_hessian_approximation_strategy()->get_apply_H()); 
                //     tr_subproblem->tr_constrained_solve(multiplication_action, g, s, box);             
                // }
                // else
                {
                    utopia_error("utopia::QuasiNewtonBound: MF solver which is not TR_box_subproblem is not suported at the moment... \n"); 
                }


                if(this->ls_strategy_) 
                    this->ls_strategy_->get_alpha(fun, g, x, s, this->alpha_);     

                s *= this->alpha_; 
                x+=s; 
                
                y = g; 
                fun.gradient(x, g);

                // norms needed for convergence check
                g_norm = this->criticality_measure_infty(x, g); 
                s_norm = norm2(s);

                // diff between fresh and old grad...
                y = g - y; 
                this->update(s, y);


                Scalar energy; 
                fun.value(x, energy); 

                // print iteration status on every iteration
                if(this->verbose_)
                    PrintInfo::print_iter_status(it, {g_norm, energy,  s_norm, this->alpha_});
                
                // check convergence and print interation info
                converged = this->check_convergence(it, g_norm, 9e9, s_norm);
                
                it++;
            }

            this->print_statistics(it); 
            return true;
        }
        
        virtual void set_parameters(const Parameters params) override
        {
            MatrixFreeNonLinearSolver<Vector>::set_parameters(params);
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

        
    private:
        Scalar alpha_;                                          /*!< Dumping parameter. */
        std::shared_ptr<LSStrategy> ls_strategy_;               /*!< Strategy used in order to obtain step \f$ \alpha_k \f$ */

    };


}
#endif //UTOPIA_QUASI_NEWTON_BOUND_HPP
