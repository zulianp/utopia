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
    class QuasiNewtonBound :    public QuasiNewtonBase<Vector>, 
                                public VariableBoundSolverInterface<Vector>

    {
        typedef UTOPIA_SCALAR(Vector)                           Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;
        
        typedef utopia::LSStrategy<Vector>                      LSStrategy;
        typedef utopia::HessianApproximation<Vector>            HessianApproximation;

        typedef utopia::MatrixFreeQPSolver<Vector>           QPSolver;
        
        
    public:
        QuasiNewtonBound(   const std::shared_ptr <HessianApproximation> &hessian_approx,
                            const std::shared_ptr <QPSolver> &linear_solver,
                            const Parameters params = Parameters()):
                            QuasiNewtonBase<Vector>(hessian_approx, linear_solver)
        {
            set_parameters(params);
        }
        
        bool solve(FunctionBase<Vector> &fun, Vector &x) override
        {
            using namespace utopia;
            
            Vector g, s, y; 
            
            Scalar g_norm, r_norm, g0_norm, s_norm=1;
            SizeType it = 0;

            Scalar alpha = 1.0; 
            
            bool converged = false;

            this->make_iterate_feasible(x); 
            
            fun.gradient(x, g);
            g0_norm = norm2(g);
            
            if(this->verbose_) {
                this->init_solver("QUASI NEWTON BOUND", {" it. ", "|| g ||", "r_norm", "|| p_k || ", "alpha"});
                PrintInfo::print_iter_status(it, {g0_norm, s_norm});
            }
            it++; 
            
            this->initialize_approximation(); 

            while(!converged)
            {
                s = local_zeros(local_size(x)); 

                if(QPSolver * qp_solver = dynamic_cast<QPSolver*>(this->linear_solver().get()))
                {
                    auto box = this->build_correction_constraints(x);
                    qp_solver->set_box_constraints(box); 
                    auto multiplication_action = this->hessian_approx_strategy_->build_apply_H(); 
                    qp_solver->solve(*multiplication_action, -1.0*g, s);             
                }
                else
                {
                    utopia_error("utopia::QuasiNewtonBound: MF solver which is not QPSolver is not suported at the moment... \n"); 
                }

                alpha = this->get_alpha(fun, g, x, s); 
                s *= alpha; 
                x+=s; 
                
                y = g; 
                fun.gradient(x, g);

                // norms needed for convergence check
                g_norm = this->criticality_measure_infty(x, g); 
                s_norm = norm2(s);
                r_norm = g_norm/g0_norm;

                // diff between fresh and old grad...
                y = g - y; 
                this->update(s, y);

                // print iteration status on every iteration
                if(this->verbose_)
                    PrintInfo::print_iter_status(it, {g_norm, r_norm,  s_norm, alpha});
                
                // check convergence and print interation info
                converged = this->check_convergence(it, g_norm, 9e9, s_norm);
                
                it++;
            }

            this->print_statistics(it); 
            return true;
        }
        
        virtual void set_parameters(const Parameters params) override
        {
            QuasiNewtonBase<Vector>::set_parameters(params);
        }
        
    };


}
#endif //UTOPIA_QUASI_NEWTON_BOUND_HPP
