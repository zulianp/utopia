#ifndef UTOPIA_SOLVER_TRUSTREGION_HPP
#define UTOPIA_SOLVER_TRUSTREGION_HPP

#include "utopia_NewtonBase.hpp"
#include "utopia_TRBase.hpp"
#include "utopia_TRSubproblem.hpp"
#include "utopia_Dogleg.hpp"
#include "utopia_SteihaugToint.hpp"


 namespace utopia
 {
    	template<class Matrix, class Vector>
     	class TrustRegion final: public NewtonBase<Matrix, Vector>,
                               public TrustRegionBase<Vector>
      {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        typedef utopia::TRSubproblem<Matrix, Vector> TRSubproblem;

        typedef utopia::TrustRegionBase<Vector> TrustRegionBase;
        typedef utopia::NewtonBase<Matrix, Vector> NonLinearSolver;

     	public:
      TrustRegion(const std::shared_ptr<TRSubproblem> &tr_subproblem): 
                  NonLinearSolver(tr_subproblem)
      {
        
      }

      using utopia::TrustRegionBase<Vector>::get_pred; 

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

      bool solve(Function<Matrix, Vector> &fun, Vector &x_k) override
      {
         using namespace utopia;

         // passing solver and parameters into subproblem
         bool converged = false;
         NumericalTollerance<Scalar> tol(this->atol(), this->rtol(), this->stol());

         Scalar delta, product, ared, pred, rho, E, E_k, E_k1;

         SizeType it = 0;
         SizeType it_successful = 0;
         
         static const auto infty = std::numeric_limits<Scalar>::infinity();
         Scalar g_norm = infty, g0_norm = infty, r_norm = infty, s_norm = infty;

         bool rad_flg = false;

         Vector g, p_k = x_k, x_k1 = x_k;
         Matrix H;


         fun.gradient(x_k, g);
          g0_norm = norm2(g);
          g_norm = g0_norm;

          #define DEBUG_mode

          // TR delta initialization
          delta =  this->delta_init(x_k , this->delta0(), rad_flg);

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
            this->init_solver("TRUST_REGION_BASE", {" it. ", "|| g ||", "r_norm", "J_k", "J_{k+1}", "rho", "delta_k", "|| p_k ||"});
            PrintInfo::print_iter_status(it, {g_norm});
          }
        #endif

        bool accepted = true; 

        // solve starts here
        while(!converged)
        {
          fun.value(x_k, E_k);
          
          if(accepted)
            fun.hessian(x_k, H);
    //----------------------------------------------------------------------------
    //     new step p_k w.r. ||p_k|| <= delta
    //----------------------------------------------------------------------------
          if(TRSubproblem * tr_subproblem = dynamic_cast<TRSubproblem*>(this->linear_solver().get()))
          {
            p_k *= 0; 
            tr_subproblem->current_radius(delta); 
            tr_subproblem->solve(H, -1.0 * g, p_k);      
            this->solution_status_.num_linear_solves++;  
          }
          else
          {
            utopia_warning("TrustRegion::Set suitable TR subproblem.... \n "); 
          }

          pred = this->get_pred(g, H, p_k);
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
            it_successful++;

          accepted = this->trial_point_acceptance(rho, E, E_k, E_k1, p_k, x_k, x_k1);

          if(accepted)
          {
            x_k = x_k1;
            fun.gradient(x_k, g);
            g_norm = norm2(g);
            r_norm = g_norm/g0_norm;

            norms2(g, p_k, g_norm, s_norm); 
            r_norm = g_norm/g0_norm;            
          }
          else
          {
            s_norm = norm2(p_k);
          }

    //----------------------------------------------------------------------------
    //    convergence check
    //----------------------------------------------------------------------------

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

        // some benchmarking
        TrustRegionBase::print_statistics(it, it_successful);      

          return true;
      }


      void set_trust_region_strategy(const std::shared_ptr<TRSubproblem> &tr_linear_solver)
      {
        NonLinearSolver::set_linear_solver(tr_linear_solver);
      }


      void set_linear_solver(const std::shared_ptr<LinearSolver<Matrix, Vector> > &linear_solver) override
      {
          utopia_error("TrustRegion:: do not use set_linear solver. Use set_trust_region_strateg... \n"); 
      }


    private:
      Scalar get_pred(const Vector & g, const Matrix & B, const Vector & p_k)
      {
        return (-1.0 * dot(g, p_k) -0.5 *dot(B * p_k, p_k));
      }

  };

}

//clean-up macros
#ifdef DEBUG_mode
#undef DEBUG_mode
#endif

#endif //UTOPIA_SOLVER_TRUSTREGION_HPP

