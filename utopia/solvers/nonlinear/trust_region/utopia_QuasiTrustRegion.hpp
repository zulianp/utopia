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
    	template<class Matrix, class Vector>
      /**
       * @brief      Base class for all TR solvers. Contains all general routines related to TR solvers.
       *             Design of class allows to provide different TR strategies in order to solve TR subproblem.
       */
     	class QuasiTrustRegion : public TrustRegion<Matrix, Vector>
      {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        typedef utopia::TRSubproblem<Matrix, Vector> TRSubproblem;
        typedef utopia::HessianApproximation<Vector>    HessianApproximation;

        using TrustRegionBase<Vector>::get_pred; 

     	public:
      QuasiTrustRegion( const std::shared_ptr<TRSubproblem> &tr_subproblem = std::make_shared<SteihaugToint<Matrix, Vector, HOMEMADE> >(),
                        const Parameters params = Parameters()) : 
                        TrustRegion<Matrix, Vector>(tr_subproblem, params), 
                        _has_hessian_approx_strategy(false)
      {
        set_parameters(params);
      }

      /* @brief      Sets the parameters.
      *
      * @param[in]  params  The parameters
      */
      virtual void set_parameters(const Parameters params) override
      {
        TrustRegion<Matrix, Vector>::set_parameters(params);
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
         bool converged = false;
         NumericalTollerance<Scalar> tol(this->atol(), this->rtol(), this->stol());

         Scalar delta, product, ared, pred, rho, E_taken, E_old, E_new;// alpha;

         SizeType it = 0;
         it_successful_ = 0;
         Scalar g_norm, g0_norm, r_norm, s_norm = std::numeric_limits<Scalar>::infinity();

         bool rad_flg = false;

         Vector g, y = local_zeros(local_size(x_k).get(0)), p_k = local_zeros(local_size(x_k).get(0)), x_trial; 

        // #define DEBUG_mode

        // TR delta initialization
        delta =  this->delta_init(x_k , this->delta0(), rad_flg);
        hessian_approx_strategy_->initialize();
            
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

        // solve starts here
        while(!converged)
        {
          fun.value(x_k, E_old);
    //----------------------------------------------------------------------------
    //     new step p_k w.r. ||p_k|| <= delta
    //----------------------------------------------------------------------------
          if(TRSubproblem * tr_subproblem = dynamic_cast<TRSubproblem*>(this->linear_solver_.get()))
          {
            auto multiplication_action = FunctionOperator<Vector>(hessian_approx_strategy_->get_apply_H()); 
            tr_subproblem->tr_constrained_solve(multiplication_action, g, p_k, delta);             
          }

          x_trial = x_k + p_k; 

          // scaling correction to fit into tr radius ... 
          s_norm = norm2(p_k);
          pred = this->get_pred(g, p_k); 

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
            it_successful_++;

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

          hessian_approx_strategy_->update(p_k, y);

    //----------------------------------------------------------------------------
    //    convergence check
    //----------------------------------------------------------------------------
          g_norm = norm2(g);
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
          _has_hessian_approx_strategy  = true;
          return true;
      }
      
      
      virtual bool has_hessian_approx()
      {
          return _has_hessian_approx_strategy;
      }
        

      virtual Scalar get_pred(const Vector & g, const Vector & p_k)
      {
        // compute tr ratio... 
        Scalar l_term = dot(g, p_k);
        Scalar qp_term = hessian_approx_strategy_->compute_uHu_dot(p_k); 
        return  (- l_term - 0.5 * qp_term); 
      }        
        

  private:
    SizeType it_successful_; 
    std::shared_ptr<HessianApproximation> hessian_approx_strategy_;
    bool _has_hessian_approx_strategy;

  };

}

#endif //UTOPIA_QUASI_TRUST_REGION_HPP

