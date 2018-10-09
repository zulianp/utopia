#ifndef UTOPIA_QUASI_NEWTON_BOUND_HPP
#define UTOPIA_QUASI_NEWTON_BOUND_HPP

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
     * @brief      The Quasi Newton solver.
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<class Matrix, class Vector>
    class QuasiNewtonBound : public QuasiNewton<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                           Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;
        
        typedef utopia::LSStrategy<Matrix, Vector>              LSStrategy;
        typedef utopia::HessianApproximation<Matrix, Vector>    HessianApproximation;
        typedef utopia::IterativeSolver<Matrix, Vector>         Solver;
        typedef utopia::BoxConstraints<Vector>                  BoxConstraints;
        
        
    public:
        QuasiNewtonBound(   const std::shared_ptr <HessianApproximation> &hessian_approx,
                            const std::shared_ptr <Solver> &linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector> >(),
                            const Parameters params = Parameters()):
        QuasiNewton<Matrix, Vector>(hessian_approx, linear_solver, params)
        {
            set_parameters(params);
        }
        
        bool solve(Function<Matrix, Vector> &fun, Vector &x) override
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
                this->init_solver("QUASI NEWTON", {" it. ", "|| g ||", "r_norm", "|| p_k || ", "alpha"});
                PrintInfo::print_iter_status(it, {g_norm, s_norm});
            }
            it++; 
            
            this->hessian_approx_strategy_->initialize(fun, x);
            
            while(!converged)
            {
                this->hessian_approx_strategy_->apply_Hinv(-1.0 * g, s); 
                
                if(this->ls_strategy_) 
                    this->ls_strategy_->get_alpha(fun, g, x, s, this->alpha_);     

                s *= this->alpha_; 
                x+=s; 
                
                y = g; 
                fun.gradient(x, g);

                // norms needed for convergence check
                g_norm = criticality_measure_infty(x, g); 
                s_norm = norm2(s);

                // diff between fresh and old grad...
                y = g - y; 
                this->hessian_approx_strategy_->update(s, y);

                // print iteration status on every iteration
                if(this->verbose_)
                    PrintInfo::print_iter_status(it, {g_norm, s_norm, this->alpha_});
                
                // check convergence and print interation info
                converged = this->check_convergence(it, g_norm, 9e9, s_norm);
                
                it++;
            }

            this->print_statistics(it); 
            return true;
        }
        
        virtual void set_parameters(const Parameters params) override
        {
            QuasiNewton<Matrix, Vector>::set_parameters(params);
        }
        
        
        virtual bool set_box_constraints(const BoxConstraints & box)
        {
            constraints_ = box;
            return true;
        }

        virtual bool  get_box_constraints(BoxConstraints & box)
        {
            box = constraints_;
            return true;
        }

    private: 
      virtual Scalar criticality_measure_infty(const Vector & x, const Vector & g)
      {

        Vector Pc; 
        Vector x_g = x - g; 
        Vector ub, lb; 

        Scalar n = local_size(x).get(0); 

        if(constraints_.has_upper_bound())
          ub = *constraints_.upper_bound(); 
        else
          ub =  local_values(n, 9e15); 

        if(constraints_.has_lower_bound())
          lb = *constraints_.lower_bound(); 
        else
          lb =  local_values(n, -9e15); 

        get_projection(x_g, lb, ub, Pc); 
        
        Pc -= x; 
        return norm2(Pc); 
      }


      bool get_projection(const Vector & x, const Vector &lb, const Vector &ub, Vector & Pc)
      {
          Pc = local_values(local_size(x).get(0), 1.0);
          {
              Read<Vector> r_ub(ub), r_lb(lb), r_x(x);
              Write<Vector> wv(Pc); 

              each_write(Pc, [ub, lb, x](const SizeType i) -> double { 
                          Scalar li =  lb.get(i); Scalar ui =  ub.get(i); Scalar xi =  x.get(i);  
                          if(li >= xi)
                            return li; 
                          else
                            return (ui <= xi) ? ui : xi; }   );
          }

          return true;
      }

    void make_iterate_feasible(Vector & x)
    {
        if(!constraints_.has_upper_bound() || !constraints_.has_lower_bound())
            return; 

        const Vector x_old = x; 

        if(constraints_.has_upper_bound() && constraints_.has_lower_bound())
        {
            const auto &ub = *constraints_.upper_bound();
            const auto &lb = *constraints_.lower_bound();

            {
              Read<Vector> r_ub(ub), r_lb(lb), r_x(x_old);
              Write<Vector> wv(x); 

              each_write(x, [ub, lb, x_old](const SizeType i) -> double { 
                          Scalar li =  lb.get(i); Scalar ui =  ub.get(i); Scalar xi =  x_old.get(i);  
                          if(li >= xi)
                            return li; 
                          else
                            return (ui <= xi) ? ui : xi; }   );
            }
        }
        else if(constraints_.has_upper_bound() && !constraints_.has_lower_bound())
        {
            const auto &ub = *constraints_.upper_bound();

            {
              Read<Vector> r_ub(ub), r_x(x_old);
              Write<Vector> wv(x); 

              each_write(x, [ub, x_old](const SizeType i) -> double { 
                          Scalar ui =  ub.get(i); Scalar xi =  x_old.get(i);  
                            return (ui <= xi) ? ui : xi; }   );
            }
        }
        else
        {
            const auto &lb = *constraints_.lower_bound();

            {
              Read<Vector> r_lb(lb), r_x(x_old);
              Write<Vector> wv(x); 

              each_write(x, [lb, x_old](const SizeType i) -> double { 
                          Scalar li =  lb.get(i); Scalar xi =  x_old.get(i);  
                          return (li >= xi) ? li : xi; }   );
            }
        }    

    }      


    
    private:
        BoxConstraints                  constraints_;

        
    };


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TODO: 
// - line solver should be delegated 
// - create common interface for constraint nonlinear solver     
// - constraints  
// criticality_measure_infty should not create inf vectors     


}
#endif //UTOPIA_QUASI_NEWTON_BOUND_HPP
