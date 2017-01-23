/*
* @Author: alenakopanicakova
* @Date:   2016-05-11
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-11-08
*/
//    .................................. WORK IN PROGRESS............................

#ifndef UTOPIA_TRUSTREGION_LOCAL_SOLVE_HPP
#define UTOPIA_TRUSTREGION_LOCAL_SOLVE_HPP
#include "utopia_GLFunction.hpp"
#include "utopia_TR_base.hpp"


namespace utopia 
{
  template<class GlobalMatrix, class GlobalVector, class LocalMatrix, class LocalVector>
  class Local_TR :  public TrustRegionBase<LocalMatrix, LocalVector> 
  {
    typedef UTOPIA_SCALAR(LocalVector)    Scalar;
    typedef UTOPIA_SIZE_TYPE(LocalVector) SizeType;
    typedef utopia::LinearSolver<LocalMatrix, LocalVector> Solver;
    typedef utopia::TRSubproblem<LocalMatrix, LocalVector> TRSubproblem; 

    public:

      Local_TR( const std::shared_ptr<TRSubproblem> &tr_subproblem = std::shared_ptr<TRSubproblem>(),
                const Parameters params = Parameters())
                                              : TrustRegionBase<LocalMatrix, LocalVector>(tr_subproblem, params)  
      {
        set_parameters(params);        
      }

      Scalar solve(GLFunction<GlobalMatrix, GlobalVector, LocalMatrix, LocalVector> &fun, 
                    const GlobalVector &x, const GlobalVector &g, LocalVector &s, 
                    const Scalar& delta_L)
      {
        using namespace utopia;
        this->delta_max(delta_L); 

        // max delta initialized from global solve 
        delta               = delta_L; 
        delta_working       = delta_L; 

        LocalVector x_k, x_0;
        LocalVector g_0, g_diff, g_k; 
        LocalMatrix H_k;
        
        SizeType it = 0;
        Scalar g_norm, r_norm, g0_norm, s_norm; 
        bool converged = false; 


        // restriction of initial quantities 
        // init needs to be here, for restriction functions to know local sizes !!! 
        fun.restrict(x, x_0);
        fun.restrict(g, g_0);
        fun.local_value(x_0, E_0); 

        g_k = g_0; 
        x_k = x_0; 

        LocalVector p_k = x_k, x_k1 = x_k;

        fun.local_hessian(x_k,  H_k); 
        fun.local_gradient(x_k, g_k);

        g_diff  = g_0 - g_k; 
        // g_k     += g_diff; 
        g_norm  = norm2(g_k); 

        #define DEBUG_mode
        // #define LS_check
        

          // print out - just to have idea how are we starting 
        #ifdef DEBUG_mode
          if(this->verbose_)
          {
             this->init_solver("TRUST_REGION_LOCAL",
                               {" it. ", "|| g ||", "r_norm", "<g, dx>", "J_k", "J_{k+1/2}", "J_{k+1}", "ared", "pred",
                                "rho", "delta_k", "delta_working",  "|| p_k || "});
            PrintInfo::print_iter_status(it, {g_norm, 9e9, 9e9, 9e9, 9e9, E_0 }); 
          }
          
        #else
          if(this->verbose_)
          {
            this->nonlinear_solve_init("TRUST_REGION_LOCAL", {" it. ", "|| g ||", "r_norm", "J_k", "J_{k+1}", "rho", "delta_k", "|| p_k ||"}); 
            PrintInfo::print_iter_status(it, {g_norm}); 
          }
        #endif


        E = E_0;  

        #ifdef LS_check
            if(TRSubproblem * tr_subproblem = dynamic_cast<TRSubproblem*>(this->linear_solver_.get()))
              tr_subproblem->current_radius(9e9);  
              this->linear_solve(H_k, g_k, p_k);

          x_k1 = x_k + p_k; 
          fun.local_value(x_k1, E); 
          fun.local_gradient(x_k1, g_k);
          // g_k += g_diff;  
          g_norm = norm2(g_k);
          it++; 
          PrintInfo::print_iter_status(it, {g_norm, 9e9, 9e9, 9e9, 9e9, E_0 }); 
        #endif


        while(!converged)
        {
          x_k = x_k1; 
          E_k = E; 


          fun.local_gradient(x_k, g_k); 
          fun.local_hessian(x_k, H_k); 
          // g_k += g_diff; 
          
      //----------------------------------------------------------------------------
      //     new step p_k w.r. ||p_k|| <= delta
      //----------------------------------------------------------------------------
          if(TRSubproblem * tr_subproblem = dynamic_cast<TRSubproblem*>(this->linear_solver_.get()))
            // tr_subproblem->current_radius(delta_working);  
            tr_subproblem->current_radius(delta);  

          this->linear_solve(H_k, g_k, p_k);
          this->get_pred(g_k, H_k, p_k, pred); 

          // value of the objective function with correction 
          fun.local_value(x_k + p_k, E_k1);
          product = dot(g_k, p_k);  // just to do tests 

          // decrease ratio 
          ared = E_k - E_k1;        // reduction observed on objective function
          rho = ared/ pred;   // decrease ratio 

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
          fun.local_gradient(x_k1, g_k);
          // g_k += g_diff;  
          g_norm = norm2(g_k);
          s_norm = norm2(p_k); 
        

          #ifdef DEBUG_mode
            if(this->verbose_)
              PrintInfo::print_iter_status(it, {g_norm, g_norm, product, E_k,  E_k1, E, ared, pred,  rho, delta, delta_working, s_norm}); 
          #else
            if(this->verbose_)
              PrintInfo::print_iter_status({it, g_norm, E_k, E, rho, delta}); 
          #endif

          converged = this->check_convergence(it, g_norm, g_norm, g_norm, delta); 
          if(E < -1e-5)
            converged = true;

    //----------------------------------------------------------------------------
    //      tr. radius update 
    //----------------------------------------------------------------------------
          this->delta_update(rho, p_k, delta);

          // if(ared > 0 && rho > this->eta2())
          //   delta_working = this->delta_max() - norm2(p_k); 

          it++; 
        }

          x_k = x_k1;                 //    - in case, we need to pass solution out - pass this one 
          s = x_k - x_0;              //    local correction 
          return E_0 - E;             //    local reduction 
      }




    void set_parameters(const Parameters params) override
    {
      TrustRegionBase<LocalMatrix, LocalVector>::set_parameters(params);
    }


  private:
    Scalar product, ared, pred, rho; 
    Scalar E_0, E, E_k, E_k1; 

    Scalar delta, delta_working; 


    virtual bool delta_update(const Scalar &rho, const LocalVector &p_k, Scalar &delta)
    {
        if(rho < this->eta1())
        {
          // in L2 norm !!! 
          // delta           = this->gamma1() * norm2(p_k); 
          delta           = this->gamma1() * delta; 
        }
        else if (rho > this->eta2())
        {
          delta         = std::min(this->gamma2() * delta, this->delta_max()); 
          // delta_working = this->delta_max() - norm2(p_k); 
        }      
          return true; 
    }


    virtual bool check_convergence(const SizeType &it, const Scalar & g_norm, const Scalar & r_norm, const Scalar & s_norm, const Scalar & delta)
    { 
      if(s_norm > this->delta_max())
      {
          std::cout<<"LocalTR:: exceeding max_l   \n"; 
          return true;       
      }

      return TrustRegionBase<LocalMatrix, LocalVector>::check_convergence(it, g_norm, r_norm, s_norm, delta); 
    }


  };

}

#endif //UTOPIA_TRUSTREGION_LOCAL_SOLVE_HPP

