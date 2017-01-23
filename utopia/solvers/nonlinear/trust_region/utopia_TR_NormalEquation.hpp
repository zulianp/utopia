/*
* @Author: alenakopanicakova
* @Date:   2016-05-18
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-10-11
*/

#ifndef UTOPIA_TRUSTREGION_NORMAL_EQ_HPP
#define UTOPIA_TRUSTREGION_NORMAL_EQ_HPP
#include "utopia_Dogleg.hpp"
#include "utopia_SteihaugToint.hpp"
#include "utopia_TR_base.hpp"
#include "utopia_TRSubproblem.hpp"
#include "utopia_CauchyPoint.hpp"
#include "utopia_NonLinearLeastSquaresSolver.hpp"

     namespace utopia 
     {
    	template<class Matrix, class Vector>
      /**
       * @brief      Trust region solver for Normal Equations. 
       *             Solution process is not very different from one in TrustRegionBase class, but it is specialized for normal eq. 
       */     
     	class LeastSquaresTrustRegion : public NonLinearLeastSquaresSolver<Matrix, Vector>,
                                      public TrustRegionBase<Matrix, Vector> 
      {
     		typedef typename utopia::Traits<Vector>::Scalar Scalar;
     		typedef utopia::LinearSolver<Matrix, Vector> Solver;
        typedef utopia::TRSubproblem<Matrix, Vector> TRSubproblem; 
        typedef utopia::NonLinearLeastSquaresSolver<Matrix, Vector> NonLinearLeastSquaresSolver;
     		typedef typename utopia::Traits<Vector>::SizeType SizeType;
        typedef TrustRegionBase<Matrix, Vector> TrustRegionBase;


     	public:
      LeastSquaresTrustRegion(const std::shared_ptr<TRSubproblem> &linear_solver = std::shared_ptr<TRSubproblem>(),
                              const Parameters params = Parameters())
                              : NonLinearLeastSquaresSolver(linear_solver, params), 
                                TrustRegionBase(params)
      {
        set_parameters(params);        
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
       * @brief      TR solution process. 
       *
       * @param      fun   The fun with nonlinear context. 
       * @param      x_k   The initial guess/ solution.
       *
       * @return     true
       */
      virtual bool solve(LeastSquaresFunction<Matrix, Vector> &fun, Vector &x_k)
      {

         using namespace utopia;
         bool converged = false; 

         NumericalTollerance<Scalar> tol(this->atol(), this->rtol(), this->stol());
         

         Scalar it = 0, rad_flg;
         Scalar g_norm, g0_norm, r_norm, s_norm = std::numeric_limits<Scalar>::infinity(); 

         Vector r_k , p_CP = x_k, dx = x_k, p_N = x_k, p_k = x_k, x_k1 = x_k, g;
         Matrix J_k, J_T, H;


         fun.jacobian(x_k, J_k); 
         fun.residual(x_k, r_k);


        #define DEBUG_mode
        //  #define LS_check


  			// TR delta initialization
        delta = this->delta0_; 
        rad_flg = this->delta_init(x_k ,delta); 


        // just to start
        // if(params().verbose()) 
        //   this->info().TR_init_message(params()); 
        g0_norm = norm2(r_k);
        g_norm = g0_norm;

        // print out - just to have idea how we are starting 
        #ifdef DEBUG_mode
          if(this->verbose_)
          {
              this->init_solver("TRUST_NORMAL_EQ",
                                {" it. ", "|| g ||", "r_norm", "<g, dx>", "J_k", "J_{k+1/2}", "J_{k+1}", "ared", "pred",
                                 "rho", "delta_k", "|| p_k || "});
            PrintInfo::print_iter_status({it, g_norm}); 
          }
          
        #else
          if(this->verbose_)
          {
            this->nonlinear_solve_init("TRUST_NORMAL_EQ", {" it. ", "|| g ||", "r_norm", "J_k", "J_{k+1}", "rho", "delta_k", "|| p_k ||"}); 
            PrintInfo::print_iter_status({it, g_norm}); 
          }
        #endif


        // found out if there is a linear solution - or start with the newton step 
        #ifdef LS_check
          if(this->linear_solution_check(fun, g, H, p_N, x_k))
          {
            if(this->verbose_) 
              TrustRegionBase::check_convergence(*this, tol, this->max_it(), 0, 0, 0, 0, 0); 
            return true; 
          }
          x_k1 = x_k; 
        #endif

        it++; 
        fun.value(x_k, E); 
        fun.residual(x_k, r_k);

        while(!converged)
        {
          E_k = E; 

          fun.jacobian(x_k, J_k); 
          J_T = transpose(J_k); 

          g = J_T * r_k; 
          H = J_T * J_k;

    //----------------------------------------------------------------------------
    //     new step p_k w.r. ||p_k|| <= delta
    //----------------------------------------------------------------------------
          // this->tr_subproblem->get_pk(r_k, g, J_k, H, delta, p_k, pred); 
        if(TRSubproblem * tr_subproblem = dynamic_cast<TRSubproblem*>(this->linear_solver_.get()))
          tr_subproblem->current_radius(delta);  

          this->linear_solve(J_k, r_k, p_k);
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
          product = dot(r_k, p_k);    // just for testing purposes 

          ared = E_k - E_k1;          // reduction observed on objective function
          rho = ared/ (2 * pred);     // decrease ratio 

    //----------------------------------------------------------------------------
    //     acceptance of trial point 
    //----------------------------------------------------------------------------
          // to avoid numerical errors 
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
          fun.residual(x_k, r_k); 
          
          g_norm = norm2(r_k); 
          r_norm = g_norm/g0_norm;
          s_norm = norm2(p_k); 

          #ifdef DEBUG_mode
            if(this->verbose_)
              PrintInfo::print_iter_status({it, g_norm, r_norm, product, E, E_k, E_k1, ared, pred, rho, delta, s_norm}); 
          #else
            if(this->verbose_)
              PrintInfo::print_iter_status({it, g_norm, r_norm, E, E_k1, rho, delta, s_norm}); 
          #endif

          converged = TrustRegionBase::check_convergence(*this, tol, this->max_it(), it, g_norm, r_norm, s_norm, delta); 

    //----------------------------------------------------------------------------
    //      tr. radius update 
    //----------------------------------------------------------------------------
          delta_update(rho, p_k, delta); 
          it++; 
        }
          return false;
      }


      /* @brief      Sets the parameters.
      *
      * @param[in]  params  The parameters
      */
      virtual void set_parameters(const Parameters params) override
      {

        NonLinearLeastSquaresSolver::set_parameters(params);
        TrustRegionBase::set_parameters(params);

        delta_max_  = params.delta_max();
        delta0_     = params.delta0();
        gamma1_     = params.gamma1();
        gamma2_     = params.gamma2();
        eta1_       = params.eta1();
        eta2_       = params.eta2();
        rho_tol_    = params.rho_tol();
        eps_        = params.eps();
      }



  private:
    Scalar delta, product, ared, pred, rho, E, E_k, E_k1; 
    std::shared_ptr<TRSubproblem> tr_subproblem;

    Scalar delta_max_; 
    Scalar delta0_; 
    Scalar gamma1_;
    Scalar gamma2_;
    Scalar eta1_;
    Scalar eta2_;
    Scalar rho_tol_;
    Scalar eps_; 
      
     /**
      * @brief      update of tr radius, specialized for solution in L2 norm. 
      *             In case u need to solve eq in || \cdot ||_{L_{\infty}}, there is need to change the 1st if statement.
      *
      * @param[in]  rho    actual/predicted reduction
      * @param      delta  new radius
      * @param      p_k    current step 
      *
      * @return    
      */
    bool delta_update(const Scalar &rho, const Vector &p_k, Scalar &delta) override
    {
        if(rho < eta1_)
        {
          Scalar dx_norm = norm2(p_k); 
          delta = gamma1_ * dx_norm;   // this is fine for L2 norm 
        }
        else if (rho > eta2_)
        {
          delta = std::min(gamma2_ * delta, delta_max_); 
        }      
        return true; 
    }

  };

}

#endif //UTOPIA_TRUSTREGION_NORMAL_EQ_HPP

