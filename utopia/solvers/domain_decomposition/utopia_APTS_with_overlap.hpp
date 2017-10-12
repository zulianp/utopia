/*
* @Author: alenakopanicakova
* @Date:   2016-05-19
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-07-03
*/
//    .................................. WORK IN PROGRESS............................

#ifndef UTOPIA_APTS_SOLVER_HPP
#define UTOPIA_APTS_SOLVER_HPP

#include "utopia_GLFunction.hpp"
#include "utopia_TRBase.hpp"
#include "utopia_TR_LocalSolve.hpp" 


#include <iomanip>
#include <iostream>
#include <limits>
#include <algorithm>


namespace utopia 
{
    template<class GlobalMatrix, class GlobalVector, class LocalMatrix, class LocalVector>
    class APTS_with_overlap :  public TrustRegionBase<GlobalMatrix, GlobalVector> 
    {
        typedef UTOPIA_SCALAR(GlobalVector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(GlobalVector) SizeType;
        typedef utopia::LinearSolver<GlobalMatrix, GlobalVector> Solver;
        typedef utopia::TRSubproblem<GlobalMatrix, GlobalVector> TRSubproblem; 


    public:
      APTS_with_overlap(  Local_TR<GlobalMatrix, GlobalVector, LocalMatrix, LocalVector> local_TR,
                          const std::shared_ptr<TRSubproblem> &tr_subproblem = std::shared_ptr<TRSubproblem>(), 
                          const Parameters params = Parameters()):     
                                                                    TrustRegionBase<GlobalMatrix, GlobalVector>(tr_subproblem, params),
                                                                    local_TR(local_TR)   
      {
        set_parameters(params);     
      }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// --------------------------------------------- APTS  ---------------------------------------------------
  ///   symbols without subscript are meant to represent global objects, 
  ///   while variables with subscript k represent objects on k-th subdomain  - change k here 
  ///   e.g.  s = \sum I_k s_k      
  ///   TODO:: change this to be compatible with other TR solvers !!!
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////


      bool solve(GLFunction<GlobalMatrix, GlobalVector, LocalMatrix, LocalVector> &fun, GlobalVector &x_k) 
      {
        using namespace utopia;

        SizeType it = 0;  
        Scalar r_norm, s_norm, g0_norm, g_norm = std::numeric_limits<Scalar>::infinity();  // norm of global gradient 

        GlobalVector g;                 // global gradient 
        GlobalVector p_N = x_k;         // Newton step on 1st iteration
        GlobalVector x_k1 = x_k;        // trial iterate 
        GlobalVector s;                 // step on given iteration
        
        GlobalMatrix H;                 // global Hessian 
        LocalVector s_k;                // local corrections
        Scalar local_red;               // reduction from local model 

        bool converged = false; 

        fun.gradient(x_k, g);
        fun.value(x_k, E); 
        g0_norm = norm2(g); 

        #define DEBUG_mode
        #define LS_check


        // print out - just to have idea how we are starting 
        #ifdef DEBUG_mode
          if(this->verbose_)
          {
              this->init_solver("APTS_no_overlap",
                                {" it. ", "|| g ||", "r_norm", "<g, dx>", "J_k", "J_{k+1/2}", "J_{k+1}", "ared", "pred",
                                 "rho", "delta_k", "|| p_k || "});
            PrintInfo::print_iter_status(it, {g0_norm, 9e9, 9e9, 9e9, 9e9, E }); 
          }
          
        #else
          if(this->verbose_)
          {
            this->init_solver("APTS_no_overlap", {" it. ", "|| g ||", "r_norm", "J_k", "J_{k+1}", "rho", "delta_k", "|| p_k ||"}); 
            PrintInfo::print_iter_status(it, {g0_norm}); 
          }
        #endif


      // found out if there is a linear solution - or start with the newton step 
      // TODO:: finish convergence check 
      #ifdef LS_check
          if(TRSubproblem * tr_subproblem = dynamic_cast<TRSubproblem*>(this->linear_solver_.get()))
            tr_subproblem->current_radius(9e9);  

          fun.hessian(x_k, H);
          this->linear_solve(H, g, p_N);
          x_k1 = x_k + p_N;

          fun.gradient(x_k1, g);
          if(norm2(g) < this->atol())
          {
            std::cout<<"APTS converged before even solving => linear solution... \n"; 
            return true; 
          }

          fun.value(x_k1, E); 
          it++; 
      #endif
        
             

        // do this more sophistificated 
        // delta = norm2(g); 
        delta = 10000; 
        delta_local = delta/ mpi_world_size();
      
      while(!converged)
      {
          x_k = x_k1; 
          E_k = E; 

          fun.gradient(x_k, g); 
          fun.hessian(x_k, H); 

          // check updates for this + max in TR_local   
          local_red = local_TR.solve(fun, x_k, g, s_k, delta_local);     // local solve 
          fun.interpolate(s_k, s);                								       // assign local corrections into glob. vector 
          s_norm = norm2(s); 

          Scalar sk_norm = norm2(s_k); 


          if(s_norm > delta)
          {
              shrinking_factor  = delta/s_norm;
              s                 *= shrinking_factor;
              delta             *= shrinking_factor;
              sk_norm           *= shrinking_factor;
              s_norm            = delta;

              fun.restrict(s, s_k);
              LocalVector local_x;
              fun.restrict(x_k, local_x);
              Scalar  v1, v2; 
              fun.local_value(local_x, v1); 
              fun.local_value(local_x + s_k, v2); 

              local_red = v1 + v2;
          }

          GlobalVector model_reduction = local_zeros(1);   
          {
            Range r = range(model_reduction);
            Write<GlobalVector> w (model_reduction);
            model_reduction.set(r.begin(), local_red);
          }
          
          // model reduction:  \sum (H_k(u_{k,0} -  H_k(u_{k,m}))
            Scalar total_model_reduction = sum(model_reduction);   


          // value of the objective function with correction 
          fun.value(x_k + s, E_k1);
          product = dot(g, s);            //  check if direction gives dicrease 

      //----------------------------------------------------------------------------
      //     acceptance of trial point 
      //----------------------------------------------------------------------------
          ared = E_k - E_k1;                              // reduction observed on objective function
          pred = std::abs(total_model_reduction);         // model reduction 
          rho = ared/ pred;                               // decrease ratio 

          if(ared < 0 || pred < 0)
            rho = 0; 
          this->trial_point_acceptance(rho, E, E_k, E_k1, s, x_k, x_k1);

      //----------------------------------------------------------------------------
      //    convergence check 
      //----------------------------------------------------------------------------
          fun.gradient(x_k1, g); 
          g_norm = norm2(g); 
          r_norm = g_norm/g0_norm;

          #ifdef DEBUG_mode
            if(this->verbose_)
              PrintInfo::print_iter_status(it, {g_norm, r_norm,  product, E_k,  E_k1, E, ared, pred,  rho, delta, s_norm});
          #else
            if(this->verbose_)
              PrintInfo::print_iter_status({it}, g_norm, r_norm, E_k, E, rho, delta, s_norm}); 
          #endif
                
           // needs to be fixed      
          converged = this->check_convergence(it, g_norm, r_norm, s_norm, delta); 

    //----------------------------------------------------------------------------
    //      tr. radius update 
    //----------------------------------------------------------------------------
          delta_update(rho, s, delta, delta_local); 
          it++; 
        }
          x_k = x_k1; 
          return true;
      }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////        

    void set_parameters(const Parameters params) override
    {
      TrustRegionBase<GlobalMatrix, GlobalVector>::set_parameters(params);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////        
    private:

      Scalar delta, delta_local, product;  
      Scalar ared, pred, rho; 
      Scalar E, E_k, E_k1; 

      Scalar shrinking_factor; 

      Local_TR<GlobalMatrix, GlobalVector, LocalMatrix, LocalVector> local_TR;       /*!< local TR   */  


    virtual bool delta_update(const Scalar &rho, const GlobalVector &p_k, Scalar &radius, Scalar & local_radius)
    {
        if(rho < this->eta1())
        {
          // in L2 norm !!! 
          radius       = this->gamma1() * norm2(p_k); 
          local_radius = this->gamma1() * local_radius;
        }
        else if (rho > this->eta2())
        {
          radius        = std::min(this->gamma2() * radius, this->delta_max()); 
          local_radius   = std::min(this->gamma2() * local_radius, this->delta_max()); 
        }      
          return true; 
    }

  };

}
#endif //UTOPIA_GASPIN_SOLVER_HPP
