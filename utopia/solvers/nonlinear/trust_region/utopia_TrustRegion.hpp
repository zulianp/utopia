/*
* @Author: alenakopanicakova
* @Date:   2016-05-11
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-11-07
*/

#ifndef UTOPIA_SOLVER_TRUSTREGION_HPP
#define UTOPIA_SOLVER_TRUSTREGION_HPP
#include "utopia_NonLinearSolver.hpp"
#include "utopia_TR_base.hpp"
#include "utopia_TRSubproblem.hpp"
#include "utopia_Dogleg.hpp"
#include "utopia_SteihaugToint.hpp"
#include "utopia_Parameters.hpp"    


 namespace utopia 
 {
    	template<class Matrix, class Vector>
      /**
       * @brief      Base class for all TR solvers. Contains all general routines related to TR solvers.
       *             Design of class allows to provide different TR strategies in order to solve TR subproblem. 
       */ 
     	class TrustRegion : public NonLinearSolver<Matrix, Vector>, 
                          public TrustRegionBase<Matrix, Vector> 
      {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
     		typedef utopia::LinearSolver<Matrix, Vector> Solver;
        typedef utopia::TRSubproblem<Matrix, Vector> TRSubproblem; 
        typedef utopia::TrustRegionBase<Matrix, Vector> TrustRegionBase; 
        typedef utopia::NonLinearSolver<Matrix, Vector> NonLinearSolver;
     	
     	public:
      TrustRegion(const std::shared_ptr<TRSubproblem> &tr_subproblem = std::shared_ptr<TRSubproblem>(),
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

      /*!
      \details
                Determine, wheater problem is linear, or we need nonlinear solve 
      @note
      \param fun          - function with evaluation routines 
      \param H            - hessian
      \param g            - gradient
      \param p_N          - Newton step
      \param x_k          - current iterate
        */
      virtual bool linear_solution_check(
        
        Function<Matrix, Vector> &fun, 
        Vector & g, 
        const Matrix & H, 
        Vector & p_N, 
        Vector & x_k)
      {
        this->linear_solve(H, -1 * g, p_N);
        fun.gradient(x_k + p_N, g);
        Scalar g_norm = norm2(g);

        if(g_norm < 1e-7)
        {
          x_k += p_N; 
          std::cout<<"To solve linear problem, TR solver is not really needed ..."; 
          return true; 
        }

        x_k += p_N; 
        return false; 
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
         Scalar rad_flg, g_norm, g0_norm, r_norm, s_norm = std::numeric_limits<Scalar>::infinity();
         Vector g, p_CP = x_k, p_N = x_k, p_k = x_k, x_k1 = x_k;
         Matrix H; 

         fun.hessian(x_k, H); 
         fun.gradient(x_k, g);

          #define DEBUG_mode
          //  #define LS_check


        // TR delta initialization
        delta = this->delta0(); 
        rad_flg = this->delta_init(x_k ,delta); 
        //delta = norm2(g);   // also possible
        // delta = 10;        // testing 

        // just to start  CHECK THIS OUT 
        // if(params().verbose()) 
        //   this->info().TR_init_message(params()); 
        
        g0_norm = norm2(g);
        g_norm = g0_norm;
        
        // print out - just to have idea how we are starting 
        #ifdef DEBUG_mode
          if(this->verbose_)
          {
              this->init_solver("TRUST_REGION_BASE",
                                {" it. ", "|| g ||", "r_norm", "<g, dx>", "J_k", "J_{k+1/2}", "J_{k+1}", "ared", "pred",
                                 "rho", "delta_k", "|| p_k || "});
            PrintInfo::print_iter_status(it, {g_norm}); 
          }
          
        #else
          if(this->verbose_)
          {
            this->nonlinear_solve_init("TRUST_REGION_BASE", {" it. ", "|| g ||", "r_norm", "J_k", "J_{k+1}", "rho", "delta_k", "|| p_k ||"}); 
            PrintInfo::print_iter_status(it, {g_norm}); 
          }
        #endif


        // found out if there is a linear solution - or start with the newton step 
        #ifdef LS_check
          if(linear_solution_check(fun, g, H, p_N, x_k))
          {
            if(this->verbose_) {
              TrustRegionBase::check_convergence(*this, tol, this->max_it(), 0, 0, 0, 0, 0); 
            }
            return true; 
          }
          x_k1 = x_k; 
        #endif

        it++; 
        fun.value(x_k, E); 
        fun.gradient(x_k, g);

        // solve starts here 
        while(!converged)
        {
          E_k = E; 
          fun.hessian(x_k, H); 
    //----------------------------------------------------------------------------
    //     new step p_k w.r. ||p_k|| <= delta
    //----------------------------------------------------------------------------          
          if(TRSubproblem * tr_subproblem = dynamic_cast<TRSubproblem*>(this->linear_solver_.get()))
            tr_subproblem->current_radius(delta);  

          this->linear_solve(H, g, p_k);
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

          this->trial_point_acceptance(rho, E, E_k, E_k1, p_k, x_k, x_k1);
    //----------------------------------------------------------------------------
    //    convergence check 
    //----------------------------------------------------------------------------
          x_k = x_k1; 
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

#endif //UTOPIA_SOLVER_TRUSTREGION_HPP

