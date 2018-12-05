#ifndef UTOPIA_QUASI_TRUST_REGION_VARIABLE_CONSTRAINT_TR_HPP
#define UTOPIA_QUASI_TRUST_REGION_VARIABLE_CONSTRAINT_TR_HPP

#include "utopia_NonLinearSolver.hpp"
#include "utopia_Parameters.hpp"    
#include "utopia_VariableBoundSolverInterface.hpp"

 namespace utopia 
 {
    	template<class Vector>
     	class QuasiTrustRegionVariableBound final:  public VariableBoundSolverInterface<Vector>, 
                                                  public TrustRegionBase<Vector>, 
                                                  public QuasiNewtonBase<Vector>
      {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        typedef utopia::MatrixFreeQPSolver<Vector>  MatrixFreeQPSolver;  
        typedef utopia::TrustRegionBase<Vector>     TrustRegionBase; 
        typedef utopia::QuasiNewtonBase<Vector>     NonLinearSolver;

        typedef utopia::HessianApproximation<Vector>      HessianApproximation;

    public:                                                                      
      QuasiTrustRegionVariableBound(const std::shared_ptr <HessianApproximation> &hessian_approx, 
                                    const std::shared_ptr<MatrixFreeQPSolver> &tr_subproblem, 
                                    const Parameters params = Parameters()) : 
                                    NonLinearSolver(hessian_approx, tr_subproblem), 
                                    it_successful_(0)

      {
        set_parameters(params);        
      }


      /* @brief      Sets the parameters.
      *
      * @param[in]  params  The parameters
      */
      void set_parameters(const Parameters params) override
      {

        NonLinearSolver::set_parameters(params);
        TrustRegionBase::set_parameters(params);
      }

      /**
       * @brief      Trust region solve. 
       *
       * @param      fun   The nonlinear solve function.
       * @param      x_k   Initial gues/ solution
       *
       *  
       * @return     true
       */
      bool solve(FunctionBase<Vector> &fun, Vector &x_k) override 
      {
        using namespace utopia;

        // passing solver and parameters into subproblem 
        bool converged = false; 
        NumericalTollerance<Scalar> tol(this->atol(), this->rtol(), this->stol());

        Scalar delta, ared, pred, rho, E_old, E_new; 

        this->make_iterate_feasible(x_k); 

        SizeType it = 0;
        it_successful_ = 0; 
        const Scalar infty = std::numeric_limits<Scalar>::infinity();
        Scalar g_norm = infty, g0_norm = infty, r_norm = infty, s_norm = infty;
        bool rad_flg = false; 

        Vector g = 0*x_k, p_k = 0*x_k, x_k1 = 0*x_k;
        Vector y =  local_zeros(local_size(x_k).get(0)); 

        fun.gradient(x_k, g);

        // TR delta initialization
        delta =  this->delta_init(x_k , this->delta0(), rad_flg); 

        this->initialize_approximation(); 
        
        g0_norm = norm2(g);
        g_norm = g0_norm;
        
        // print out - just to have idea how we are starting 
        if(this->verbose_)
        {
          this->init_solver("TRUST_REGION_BASE",
                              {" it. ", "||P_c(x-g)-x||","J_k", "J_{k+1}", "ared","pred", "rho", "delta_k", "|| p_k || "});
          PrintInfo::print_iter_status(it, {g_norm}); 
        }

        it++; 
        fun.value(x_k, E_old); 
        fun.gradient(x_k, g);


        auto multiplication_action = this->hessian_approx_strategy_->build_apply_H(); 

        // solve starts here 
        while(!converged)
        {

    //----------------------------------------------------------------------------
    //     new step p_k w.r. ||p_k|| <= delta
    //----------------------------------------------------------------------------          
          if(MatrixFreeQPSolver * tr_subproblem = dynamic_cast<MatrixFreeQPSolver*>(this->linear_solver().get()))
          {
            p_k = 0 * p_k; 
            auto box = this->merge_pointwise_constraints_with_uniform_bounds(x_k, -1.0 * delta, delta);
            tr_subproblem->set_box_constraints(box); 
            tr_subproblem->solve(*multiplication_action, -1.0*g, p_k);      
          }
          else
          {
            utopia_warning("QUasiTrustRegionVariableBound::Set suitable TR subproblem.... \n "); 
          }

          pred = this->get_pred(g, *multiplication_action, p_k); 
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
          // trial point 
          x_k1 = x_k + p_k; 

          // value of the objective function with correction 
          fun.value(x_k1, E_new);
          fun.value(x_k, E_old);

          // decrease ratio 
          ared = E_old - E_new;           // reduction observed on objective function
          pred = std::abs(pred); 
          rho = ared/ pred;               // decrease ratio         

    //----------------------------------------------------------------------------
    //     acceptance of trial point
    //----------------------------------------------------------------------------
          if(ared < 0 || pred < 0)
          {
            rho = 0;
          }
          else if(ared == pred)
          {
            rho = 1;
          }

          if (rho >= this->rho_tol())
            it_successful_++;

          // good reduction, accept trial point 
          if (rho >= this->rho_tol())
          {
            x_k += p_k;      
            
            y = g; 
            fun.gradient(x_k, g);
            y = g - y;       
          }
          // otherwise, keep old point
          else
          {
            
            Vector grad_trial; 
            fun.gradient(x_k1, grad_trial);
            y = grad_trial - g;       
          }

          this->update(p_k, y);
    //----------------------------------------------------------------------------
    //    convergence check 
    //----------------------------------------------------------------------------          
          g_norm = this->criticality_measure_infty(x_k, g); 
          s_norm = norm2(p_k); 

          if(this->verbose_)
            PrintInfo::print_iter_status(it, {g_norm, E_old, E_new, ared, pred, rho, delta, s_norm}); 

          converged = TrustRegionBase::check_convergence(*this, tol, this->max_it(), it, g_norm, r_norm, s_norm, delta); 

    //----------------------------------------------------------------------------
    //      tr. radius update 
    //----------------------------------------------------------------------------
          this->delta_update(rho, p_k, delta, true); 

          it++; 
        }
          return false;
      }


    void set_trust_region_strategy(const std::shared_ptr<MatrixFreeQPSolver> &tr_linear_solver)
    {
      NonLinearSolver::set_linear_solver(tr_linear_solver); 
    }

    private:
      SizeType it_successful_; 

  };

}

#endif //UTOPIA_QUASI_TRUST_REGION_VARIABLE_CONSTRAINT_TR_HPP

