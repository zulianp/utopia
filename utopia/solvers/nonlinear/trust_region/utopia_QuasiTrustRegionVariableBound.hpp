#ifndef UTOPIA_QUASI_TRUST_REGION_VARIABLE_CONSTRAINT_TR_HPP
#define UTOPIA_QUASI_TRUST_REGION_VARIABLE_CONSTRAINT_TR_HPP

#include "utopia_NonLinearSolver.hpp"
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
                                    const std::shared_ptr <MatrixFreeQPSolver> &tr_subproblem) : 
                                    NonLinearSolver(hessian_approx, tr_subproblem), 
                                    it_successful_(0), 
                                    initialized_(false), 
                                    loc_size_(0)

      {
        
      }

      void read(Input &in) override
      {
        TrustRegionBase::read(in); 
        NonLinearSolver::read(in); 
      }

      void print_usage(std::ostream &os) const override
      {
        TrustRegionBase::print_usage(os); 
        NonLinearSolver::print_usage(os); 
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

        this->fill_empty_bounds(local_size(x_k)); 
        this->make_iterate_feasible(x_k); 

        SizeType it = 0;
        it_successful_ = 0; 
        const Scalar infty = std::numeric_limits<Scalar>::infinity();
        Scalar g_norm = infty, g0_norm = infty, r_norm = infty, s_norm = infty;
        bool rad_flg = false; 


        SizeType loc_size_x = local_size(x_k); 
        if(!initialized_ || !g.comm().conjunction(loc_size_ == loc_size_x)) 
        {
            init_vectors(loc_size_x);
        }

        fun.gradient(x_k, g);

        // TR delta initialization
        delta =  this->delta_init(x_k , this->delta0(), rad_flg); 

        this->initialize_approximation(x_k, g); 
        auto multiplication_action = this->hessian_approx_strategy_->build_apply_H(); 


        // Vector help;
        // multiplication_action->apply(g, help);
        // delta = (norm2(g)*norm2(g))/dot(g, help); 


        g0_norm = this->criticality_measure_infty(x_k, g); 
        g_norm = g0_norm;
        
        // print out - just to have idea how we are starting 
        if(this->verbose_)
        {
          this->init_solver("QUasiTrustRegionVariableBound",
                              {" it. ", "||P_c(x-g)-x||","J_k", "J_{k+1}", "ared","pred", "rho", "delta_k", "|| p_k ||_{inf} "});
          PrintInfo::print_iter_status(it, {g_norm}); 
        }

        it++; 
        fun.value(x_k, E_old); 


        // solve starts here 
        while(!converged)
        {
          // ----------------------------------------------------------------------------
          //     new step p_k w.r. ||p_k|| <= delta
          // ----------------------------------------------------------------------------          
          if(MatrixFreeQPSolver * tr_subproblem = dynamic_cast<MatrixFreeQPSolver*>(this->linear_solver().get()))
          {
            p_k.set(0.0); 
            auto box = this->merge_pointwise_constraints_with_uniform_bounds(x_k, -1.0 * delta, delta);
            tr_subproblem->set_box_constraints(box); 
            g_help = -1.0*g; 
            tr_subproblem->solve(*multiplication_action, g_help, p_k);     
            
            std::cout<<"delta: "<< delta << " \n";

            std::cout<<"\n norm2(p_k) :" << norm2(p_k) << "  \n"; 
            std::cout<<"\n norm_infty(p_k) :" << norm_infty(p_k) << "  \n"; 

            this->solution_status_.num_linear_solves++;  
          }
          else
          {
            utopia_warning("QUasiTrustRegionVariableBound::Set suitable TR subproblem.... \n "); 
          }

          pred = this->get_pred(g, *multiplication_action, p_k);    
          // ----------------------------------------------------------------------------
          // ----------------------------------------------------------------------------
          // trial point 
          x_trial = x_k + p_k; 
          

          // value of the objective function with correction 
          fun.value(x_trial, E_new);

          // decrease ratio 
          ared = E_old - E_new;           // reduction observed on objective function
          pred = std::abs(pred); 
          rho = ared/ pred;               // decrease ratio         
  
    //----------------------------------------------------------------------------
    //     acceptance of trial point
    //----------------------------------------------------------------------------
          if(ared < 0.0 || pred < 0.0)
          {
            rho = 0.0;
          }
          else if(ared == pred)
          {
            rho = 1.0;
          }

          if (rho >= this->rho_tol()){
            it_successful_++;
          }


          // good reduction, accept trial point 
          if (rho >= this->rho_tol())
          {
            x_k = x_trial;   
            E_old = E_new; 
            
            y = g; 
            fun.gradient(x_k, g);
            y = g - y;       
            
          }
          // otherwise, keep old point
          else
          {
            fun.gradient(x_trial, g_help);
            y = g_help - g;   
          }

          this->update(p_k, y, x_k, g);
    //----------------------------------------------------------------------------
    //    convergence check 
    //----------------------------------------------------------------------------          
          g_norm = this->criticality_measure_infty(x_k, g); 
          s_norm = norm_infty(p_k); 

          if(this->verbose_){
            PrintInfo::print_iter_status(it, {g_norm, E_old, E_new, ared, pred, rho, delta, s_norm}); 
          }

          converged = TrustRegionBase::check_convergence(*this, tol, this->max_it(), it, g_norm, r_norm, s_norm, delta); 
    //----------------------------------------------------------------------------
    //      tr. radius update 
    //----------------------------------------------------------------------------
          this->delta_update_inf(rho, p_k, delta); 

          it++; 
        }
          return false;
      }


    void set_trust_region_strategy(const std::shared_ptr<MatrixFreeQPSolver> &tr_linear_solver)
    {
      NonLinearSolver::set_linear_solver(tr_linear_solver); 
    }


    private:
        void init_vectors(const SizeType &ls)
        {
            auto zero_expr = local_zeros(ls);

            p_k     = zero_expr;
            g       = zero_expr;
            g_help = zero_expr;
            y       = zero_expr;
            x_trial = zero_expr;
                           
            initialized_ = true;    
            loc_size_ = ls;                                        
        }


      SizeType it_successful_; 
      Vector g, y, p_k, x_trial, g_help;
      bool initialized_; 
      SizeType loc_size_;         

  };

}

#endif //UTOPIA_QUASI_TRUST_REGION_VARIABLE_CONSTRAINT_TR_HPP

