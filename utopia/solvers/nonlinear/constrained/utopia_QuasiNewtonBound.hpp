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
                            const std::shared_ptr <QPSolver> &linear_solver):
                            QuasiNewtonBase<Vector>(hessian_approx, linear_solver), 
                            initialized_(false), 
                            loc_size_(0)
        {

        }

        bool solve(FunctionBase<Vector> &fun, Vector &x) override
        {
            using namespace utopia;

            Scalar g_norm, r_norm, g0_norm, s_norm=1;
            SizeType it = 0;

            Scalar alpha = 1.0;
            bool converged = false;

            SizeType loc_size_x = local_size(x); 

            this->fill_empty_bounds(loc_size_x); 
            this->make_iterate_feasible(x);

            if(!initialized_ || !x.comm().conjunction(loc_size_ == loc_size_x)) 
            {
                init_memory(loc_size_x);
            }            

            fun.gradient(x, g);
            g0_norm = this->criticality_measure_infty(x, g);

            QuasiNewtonBase<Vector>::init_memory(x, g); 

            if(this->verbose_) {
                this->init_solver("QUASI NEWTON BOUND", {" it. ", "|| g ||", "r_norm", "|| p_k || ", "alpha"});
                PrintInfo::print_iter_status(it, {g0_norm, s_norm});
            }
            it++;

            this->initialize_approximation(x, g);
            auto multiplication_action = this->hessian_approx_strategy_->build_apply_H();

            UTOPIA_NO_ALLOC_BEGIN("Quasi_NewtonBound");
            while(!converged)
            {
                if(QPSolver * qp_solver = dynamic_cast<QPSolver*>(this->linear_solver().get()))
                {
                    auto box = this->build_correction_constraints(x);
                    qp_solver->set_box_constraints(box);
                    s.set(0.0);
                    g_minus = -1.0 * g; 
                    qp_solver->solve(*multiplication_action, g_minus, s);
                }
                else
                {
                    utopia_error("utopia::QuasiNewtonBound: MF solver which is not QPSolver is not suported at the moment... \n");
                }

                // UTOPIA_NO_ALLOC_BEGIN("Quasi Newton Bound :3");
                alpha = this->get_alpha(fun, g, x, s);
                s *= alpha;

                x = x + s;
                this->make_iterate_feasible(x);

                y = g;
                // UTOPIA_NO_ALLOC_END();


                fun.gradient(x, g);

                // UTOPIA_NO_ALLOC_BEGIN("Quasi Newton Bound :4");
                // norms needed for convergence check
                g_norm = this->criticality_measure_infty(x, g);
                s_norm = norm_infty(s);
                r_norm = g_norm/g0_norm;

                // diff between fresh and old grad...
                y = g - y;
                this->update(s, y, x, g);
                // UTOPIA_NO_ALLOC_END();

                // print iteration status on every iteration
                if(this->verbose_)
                    PrintInfo::print_iter_status(it, {g_norm, r_norm,  s_norm, alpha});

                // check convergence and print interation info
                converged = this->check_convergence(it, g_norm, 9e9, s_norm);

                it++;
            }
            UTOPIA_NO_ALLOC_END();

            this->print_statistics(it);
            return true;
        }


    private:
        void init_memory(const SizeType &ls)
        {
            auto zero_expr = local_zeros(ls);

            s       = zero_expr;
            g       = zero_expr;
            y       = zero_expr;
            g_minus = zero_expr;

            VariableBoundSolverInterface<Vector>::init_memory(ls); 
                           
            initialized_ = true;    
            loc_size_ = ls;                                        
        }


        Vector g, s, y, g_minus;
        bool initialized_; 
        SizeType loc_size_;            

    };


}
#endif //UTOPIA_QUASI_NEWTON_BOUND_HPP
