#ifndef UTOPIA_QUASI_TRUST_REGION_VARIABLE_CONSTRAINT_TR_HPP
#define UTOPIA_QUASI_TRUST_REGION_VARIABLE_CONSTRAINT_TR_HPP

#include "utopia_NonLinearSolver.hpp"
#include "utopia_TRBoxSubproblem.hpp"
#include "utopia_Parameters.hpp"    
#include "utopia_VariableBoundSolverInterface.hpp"

 namespace utopia 
 {
    	template<class Matrix, class Vector>
      /**
       * @brief      Trust region solver taking into account also bound constraints.
       */ 
     	class QuasiTrustRegionVariableBound :   public VariableBoundSolverInterface<Matrix, Vector>, 
                                              public TrustRegionBase<Matrix, Vector>, 
                                              public NonLinearSolver<Matrix, Vector>
      {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        typedef utopia::TRBoxSubproblem<Matrix, Vector>       TRBoxSubproblem;  
        typedef utopia::TrustRegionBase<Matrix, Vector>       TrustRegionBase; 
        typedef utopia::NonLinearSolver<Matrix, Vector>       NonLinearSolver;
        typedef utopia::LinearSolver<Matrix, Vector>                Solver;

        typedef utopia::HessianApproximation<Matrix, Vector>    HessianApproximation;
        typedef utopia::MatrixFreeSolverInterface<Matrix, Vector>   MFInterface;
     	
     	public:                                                                      
      QuasiTrustRegionVariableBound(const std::shared_ptr <HessianApproximation> &hessian_approx, 
                                    const std::shared_ptr<TRBoxSubproblem> &tr_subproblem, 
                                    const Parameters params = Parameters()) : 
                                    NonLinearSolver(tr_subproblem, params), 
                                    it_successful_(0),
                                    hessian_approx_strategy_(hessian_approx) 

      {

        if(Solver * mf_solver = dynamic_cast<Solver*>(this->linear_solver_.get()))
            std::cout<<"mf solver--- \n"; 
        else
          utopia_error("QuasiTrustRegionVariableBound, linear solver is missing MatrixFreeSolverInterface\n"); 

        set_parameters(params);        
      }

      /* @brief      Sets the parameters.
      *
      * @param[in]  params  The parameters
      */
      virtual void set_parameters(const Parameters params) override
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
      bool solve(Function<Matrix, Vector> &fun, Vector &x_k) override 
      {
        using namespace utopia;

        if(!hessian_approx_strategy_)
          utopia_error("You need to set hessian approx strategy before runing Quasi newton type of solver... \n"); 

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

        hessian_approx_strategy_->initialize(fun, x_k);
        
        if(MFInterface * mf_solver = dynamic_cast<MFInterface*>(this->linear_solver_.get()))
            mf_solver->initialize(hessian_approx_strategy_); 

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

        // solve starts here 
        while(!converged)
        {

    //----------------------------------------------------------------------------
    //     new step p_k w.r. ||p_k|| <= delta
    //----------------------------------------------------------------------------          
          p_k = 0 * p_k; 
          
          if(TRBoxSubproblem * tr_subproblem = dynamic_cast<TRBoxSubproblem*>(this->linear_solver_.get()))
          {
            p_k = 0 * p_k; 
            auto box = this->merge_pointwise_constraints_with_uniform_bounds(x_k, -1.0 * delta, delta); 
            tr_subproblem->tr_constrained_solve(g, p_k, box);
          }

          Scalar l_term = dot(g, p_k);

          Scalar qp_term  = 0.0; 
          if(MFInterface * mf_solver = dynamic_cast<MFInterface*>(this->linear_solver_.get()))
            qp_term = mf_solver->compute_uHu_dot(p_k); 


          // Scalar qp_term = hessian_approx_strategy_->compute_uHu_dot(p_k); 
          pred = - l_term - 0.5 * qp_term; 
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
          // trial point 
          x_k1 = x_k + p_k; 

          // value of the objective function with correction 
          fun.value(x_k1, E_new);
          fun.value(x_k, E_old);

          // decrease ratio 
          ared = E_old - E_new;                // reduction observed on objective function
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


          hessian_approx_strategy_->update(p_k, y);
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

        // some benchmarking 
        this->print_statistics(it);      

          return false;
      }


      virtual void set_linear_solver(const std::shared_ptr<LinearSolver<Matrix, Vector> > &ls) override
      {
          auto linear_solver = this->linear_solver(); 
          if (dynamic_cast<TRBoxSubproblem *>(linear_solver.get()) != nullptr)
          {
              TRBoxSubproblem * tr_sub = dynamic_cast<TRBoxSubproblem *>(linear_solver.get());
              tr_sub->set_linear_solver(ls);
          }
      }

      virtual void set_trust_region_strategy(const std::shared_ptr<TRBoxSubproblem> &tr_linear_solver)
      {
        NonLinearSolver::set_linear_solver(tr_linear_solver); 
      }


      /**
       * @brief      Sets strategy for computing step-size.
       *
       * @param[in]  strategy  The line-search strategy.
       *
       * @return
       */
      virtual void set_hessian_approximation_strategy(const std::shared_ptr<HessianApproximation> &strategy)
      {
          hessian_approx_strategy_      = strategy;
      }
      


    private:
      SizeType it_successful_; 
      std::shared_ptr<HessianApproximation> hessian_approx_strategy_;


  };

}

#endif //UTOPIA_QUASI_TRUST_REGION_VARIABLE_CONSTRAINT_TR_HPP

