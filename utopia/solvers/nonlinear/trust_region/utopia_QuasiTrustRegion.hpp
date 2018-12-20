#ifndef UTOPIA_QUASI_TRUST_REGION_HPP
#define UTOPIA_QUASI_TRUST_REGION_HPP

#include "utopia_NonLinearSolver.hpp"
#include "utopia_TRBase.hpp"
#include "utopia_TRSubproblem.hpp"
#include "utopia_Dogleg.hpp"
#include "utopia_SteihaugToint.hpp"
#include "utopia_Parameters.hpp"
#include "utopia_SteihaugToint.hpp"


 namespace utopia
 {
    	template<class Vector>
     	class QuasiTrustRegion final: public TrustRegionBase<Vector>, public QuasiNewtonBase<Vector> 
      {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        typedef utopia::MatrixFreeTRSubproblem<Vector> TRSubproblem;
        typedef utopia::HessianApproximation<Vector> HessianApproximation;
        typedef utopia::QuasiNewtonBase<Vector>     NonLinearSolver;

     	public:
      QuasiTrustRegion( const std::shared_ptr <HessianApproximation> &hessian_approx, 
                        const std::shared_ptr<TRSubproblem> &tr_subproblem): 
                        NonLinearSolver(hessian_approx, tr_subproblem)
      {
        
      }

      void read(Input &in) override
      {
        TrustRegionBase<Vector>::read(in); 
        QuasiNewtonBase<Vector>::read(in); 
      }

      void print_usage(std::ostream &os) const override
      {
        TrustRegionBase<Vector>::print_usage(os); 
        QuasiNewtonBase<Vector>::print_usage(os); 
      }



      /**
       * @brief      QUasi trust region solve.
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

         Scalar delta, product, ared, pred, rho, E_taken, E_old, E_new;// alpha;

         SizeType it = 0;
         SizeType it_successful = 0;
         Scalar g_norm, g0_norm, r_norm, s_norm = std::numeric_limits<Scalar>::infinity();

         bool rad_flg = false;

         Vector g, y = local_zeros(local_size(x_k).get(0)), p_k = local_zeros(local_size(x_k).get(0)), x_trial; 

        // #define DEBUG_mode

        // TR delta initialization
        delta =  this->delta_init(x_k , this->delta0(), rad_flg);
        this->initialize_approximation(); 
            
        fun.gradient(x_k, g);
        g0_norm = norm2(g);
        g_norm = g0_norm;


        // print out - just to have idea how we are starting
        #ifdef DEBUG_mode
          if(this->verbose_)
          {
              this->init_solver("QUASI TRUST REGION",
                                {" it. ", "|| g ||", "r_norm", "<g, dx>", "J_k", "J_{k+1/2}", "J_{k+1}", "ared", "pred",
                                 "rho", "delta_k", "|| p_k || "});
            PrintInfo::print_iter_status(it, {g_norm});
          }

        #else
          if(this->verbose_)
          {
            this->init_solver("QUASI TRUST REGION", {" it. ", "|| g ||", "r_norm", "J_k", "J_{k+1}", "rho", "delta_k", "|| p_k ||"});
            PrintInfo::print_iter_status(it, {g_norm});
          }
        #endif

        it++;
        auto multiplication_action = this->hessian_approx_strategy_->build_apply_H(); 

        // solve starts here
        while(!converged)
        {
          fun.value(x_k, E_old);
    //----------------------------------------------------------------------------
    //     new step p_k w.r. ||p_k|| <= delta
    //----------------------------------------------------------------------------
          if(TRSubproblem * tr_subproblem = dynamic_cast<TRSubproblem*>(this->linear_solver().get()))
          {
            p_k *= 0; 
            tr_subproblem->current_radius(delta); 
            tr_subproblem->solve(*multiplication_action, -1.0 * g, p_k);   
          }
          else
          {
            utopia_warning("TrustRegion::Set suitable TR subproblem.... \n "); 
          }


          x_trial = x_k + p_k; 
          pred = this->get_pred(g, *multiplication_action, p_k); 

    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
          if(it == 1 && rad_flg)
          {
            delta = norm2(p_k);
            delta *= 0.2;
          }

          // value of the objective function with correction
          fun.value(x_trial, E_new);
          product = dot(g, p_k);            // just to do tests

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
            it_successful++;

          // good reduction, accept trial point 
          if (rho >= this->rho_tol())
          {
            x_k += p_k;

            y = g; 
            fun.gradient(x_k, g);
            y = g - y; 

            E_taken = E_new; 
            
          }
          // otherwise, keep old point
          else
          {
            Vector grad_trial; 
            fun.gradient(x_trial, grad_trial);
            y = grad_trial - g;       

            E_taken = E_old;
          }

          this->update(p_k, y);

    //----------------------------------------------------------------------------
    //    convergence check
    //----------------------------------------------------------------------------
          norms2(g, p_k, g_norm, s_norm); 
          r_norm = g_norm/g0_norm;                

          #ifdef DEBUG_mode
            if(this->verbose_)
              PrintInfo::print_iter_status(it, {g_norm, r_norm, product, E_taken, E_old, E_new, ared, pred, rho, delta, s_norm});
          #else
            if(this->verbose_)
              PrintInfo::print_iter_status(it, {g_norm, r_norm, E_taken, E_new, rho, delta, s_norm});
          #endif

            converged = TrustRegionBase<Vector>::check_convergence(*this, tol, this->max_it(), it, g_norm, r_norm, 9e9, delta);
    //----------------------------------------------------------------------------
    //      tr. radius update
    //----------------------------------------------------------------------------
          this->delta_update(rho, p_k, delta);
          it++;

        }
        
          return true;
      }


    void set_trust_region_strategy(const std::shared_ptr<TRSubproblem> &tr_linear_solver)
    {
      NonLinearSolver::set_linear_solver(tr_linear_solver); 
    }


  };

}

#endif //UTOPIA_QUASI_TRUST_REGION_HPP

