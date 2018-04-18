/*
* @Author: alenakopanicakova
* @Date:   2016-05-11
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-07-03
*/

#ifndef UTOPIA_SOLVER_BOX_CONSTRAINT_TR_HPP
#define UTOPIA_SOLVER_BOX_CONSTRAINT_TR_HPP
#include "utopia_NonLinearSolver.hpp"
#include "utopia_TRBoxBase.hpp"
#include "utopia_TRBoxSubproblem.hpp"
#include "utopia_Parameters.hpp"    


 namespace utopia 
 {
    	template<class Matrix, class Vector>
      /**
       * @brief      Trust region solver taking into account also bound constrains.
       */ 
     	class TrustRegionVariableBound :  public NonLinearSolver<Matrix, Vector>, 
                                        public TrustRegionBoxBase<Matrix, Vector> 
      {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        typedef utopia::TRBoxSubproblem<Matrix, Vector> TRBoxSubproblem; 
        typedef utopia::TrustRegionBoxBase<Matrix, Vector> TrustRegionBase; 
        typedef utopia::NonLinearSolver<Matrix, Vector> NonLinearSolver;
     	
     	public:
      TrustRegionVariableBound(const std::shared_ptr<TRBoxSubproblem> &tr_subproblem = std::shared_ptr<TRBoxSubproblem>(),
                              const Parameters params = Parameters())
                              : NonLinearSolver(tr_subproblem, params)  
      {

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

        // passing solver and parameters into subproblem 
        bool converged = false, accepted = true; 
        NumericalTollerance<Scalar> tol(this->atol(), this->rtol(), this->stol());

        Scalar delta, ared, pred, rho, E_old, E_new; 

        SizeType it = 0; 
        Scalar g_norm, g0_norm, r_norm, s_norm = std::numeric_limits<Scalar>::infinity();
        bool rad_flg = false; 

        Vector g = 0*x_k, p_k = 0*x_k, x_k1 = 0*x_k;
        Matrix H; 

        fun.gradient(x_k, g);

        // TR delta initialization
        delta =  this->delta_init(x_k , this->delta0(), rad_flg); 
        // delta = 10;        // testing 

        g0_norm = norm2(g);
        g_norm = g0_norm;
        
        // print out - just to have idea how we are starting 
        if(this->verbose_)
        {
          this->init_solver("TRUST_REGION_BASE",
                              {" it. ", "||P_c(x-g)-x||","J_k", "J_{k+1}", "rho", "ared","pred", "delta_k", "|| p_k || "});
          PrintInfo::print_iter_status(it, {g_norm}); 
        }

        it++; 
        fun.value(x_k, E_old); 
        fun.gradient(x_k, g);

        // solve starts here 
        while(!converged)
        {
          if(accepted)
          {
            fun.value(x_k, E_old); 
            fun.hessian(x_k, H); 
          }
    //----------------------------------------------------------------------------
    //     new step p_k w.r. ||p_k|| <= delta
    //----------------------------------------------------------------------------          
          if(TRBoxSubproblem * tr_subproblem = dynamic_cast<TRBoxSubproblem*>(this->linear_solver_.get()))
          {
            p_k = 0 * p_k; 
            tr_subproblem->current_radius(delta);  

            Vector ub, lb; 
            this->merge_tr_with_pointwise_constrains(x_k, delta, ub, lb); 
            
            auto box = make_box_constaints(make_ref(lb), make_ref(ub)); 
            tr_subproblem->tr_constrained_solve(H, g, p_k, box);
          }

          this->get_pred(g, H, p_k, pred); 
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
          // trial point 
          x_k1 = x_k + p_k; 

          // value of the objective function with correction 
          fun.value(x_k1, E_new);

          // decrease ratio 
          ared = E_old - E_new;                // reduction observed on objective function
          pred = std::abs(pred); 
          rho = ared/ pred;               // decrease ratio         

    //----------------------------------------------------------------------------
    //     acceptance of trial point 
    //----------------------------------------------------------------------------
          if(ared < 0 || pred < 0)
            rho = 0.0; 
          else if(rho != rho)
            rho = 0.0; 


          this->trial_point_acceptance(rho, x_k1, x_k); 
    //----------------------------------------------------------------------------
    //    convergence check 
    //----------------------------------------------------------------------------
          if(accepted)
            fun.gradient(x_k, g); 
          
          g_norm = this->criticality_measure_infty(x_k, g); 
          s_norm = norm2(p_k); 

          if(this->verbose_)
            PrintInfo::print_iter_status(it, {g_norm, E_old, E_new, rho, ared, pred, delta, s_norm}); 

          converged = TrustRegionBase::check_convergence(*this, tol, this->max_it(), it, g_norm, r_norm, s_norm, delta); 

    //----------------------------------------------------------------------------
    //      tr. radius update 
    //----------------------------------------------------------------------------
          this->delta_update(rho, p_k, delta); 
          it++; 
        }

          return false;
      }



  };

}

#endif //UTOPIA_SOLVER_BOX_CONSTRAINT_TR_HPP

