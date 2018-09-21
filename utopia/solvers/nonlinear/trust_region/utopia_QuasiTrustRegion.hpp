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
        typedef utopia::HessianApproximation<Matrix, Vector>    HessianApproximation;

     	public:
      QuasiTrustRegion( const std::shared_ptr<TRSubproblem> &tr_subproblem = std::make_shared<SteihaugToint<Matrix, Vector>>(),
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

         Scalar delta, product, ared, pred, rho, E, E_k, E_k1;

         SizeType it = 0;
         it_successful_ = 0;
         Scalar g_norm, g0_norm, r_norm, s_norm = std::numeric_limits<Scalar>::infinity();

         bool rad_flg = false;

         Vector g; 
         Vector p_k = local_zeros(local_size(x_k).get(0)); 
         Matrix H;

        #define DEBUG_mode

        // TR delta initialization
        delta =  this->delta_init(x_k , this->delta0(), rad_flg);


        if(has_hessian_approx())
            hessian_approx_strategy_->initialize(fun, x_k, H);
        else
        {
            std::cout<<"hessian approx is missing..... \n"; 
            fun.hessian(x_k, H);
        }
            

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
          fun.value(x_k, E_k);
    //----------------------------------------------------------------------------
    //     new step p_k w.r. ||p_k|| <= delta
    //----------------------------------------------------------------------------
          if(TRSubproblem * tr_subproblem = dynamic_cast<TRSubproblem*>(this->linear_solver_.get()))
            tr_subproblem->tr_constrained_solve(H, g, p_k, delta);

          this->get_pred(g, H, p_k, pred);
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
          if(it == 1 && rad_flg)
          {
            delta = norm2(p_k);
            delta *= 0.2;
          }

          // value of the objective function with correction
          fun.value(x_k + p_k, E_k1);
          product = dot(g, p_k);            // just to do tests

          // decrease ratio
          ared = E_k - E_k1;                // reduction observed on objective function
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
            E = E_k1; 
          }
          // otherwise, keep old point
          else
          {
            E = E_k; 
          }

          if(has_hessian_approx())
              hessian_approx_strategy_->approximate_hessian(fun, x_k, p_k, H,  g);
          else
          {
            std::cout<<"hessian approx is missing..... \n"; 
              fun.hessian(x_k, H);
          }  

    //----------------------------------------------------------------------------
    //    convergence check
    //----------------------------------------------------------------------------
          fun.gradient(x_k, g);
          g_norm = norm2(g);
          r_norm = g_norm/g0_norm;
          s_norm = norm2(p_k);

          #ifdef DEBUG_mode
            if(this->verbose_)
              PrintInfo::print_iter_status(it, {g_norm, r_norm, product, E, E_k, E_k1, ared, pred, rho, delta, s_norm});
          #else
            if(this->verbose_)
              PrintInfo::print_iter_status(it, {g_norm, r_norm, E, E_k1, rho, delta, s_norm});
          #endif

            converged = TrustRegionBase<Matrix, Vector>::check_convergence(*this, tol, this->max_it(), it, g_norm, r_norm, 9e9, delta);
    //----------------------------------------------------------------------------
    //      tr. radius update
    //----------------------------------------------------------------------------
          this->delta_update(rho, p_k, delta);
          it++;

        }

        // some benchmarking
        this->print_statistics(it);

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
        
        

  private:
    SizeType it_successful_; 
    std::shared_ptr<HessianApproximation> hessian_approx_strategy_;
    bool _has_hessian_approx_strategy;

  };

}

#endif //UTOPIA_QUASI_TRUST_REGION_HPP

